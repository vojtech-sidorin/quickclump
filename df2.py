#!/usr/bin/env python
# -*- coding: utf-8 -*-
# script by Vojtech Sidorin

"""
Identifies clumps within a 3D FITS datacube.

The algorithm was originally conceived by Richard Wünsch, who also
published its first implementation in Python, later rewritten in C.
Compared to the original, this implementation (df2) uses different
data structures and doesn't use the parameter Nlevels.  This new
implementation is also significantly less computationally expensive
and scales linearly with the datacube volume (number of pixels).

Type "python df2.py -h" for usage help.

See http://galaxy.asu.cas.cz/~richard/dendrofind/ for the description
of the original algorithm.  First practical usage together with another
description was published by Wünsch et al. (2012)
(http://adsabs.harvard.edu/abs/2012A%26A...539A.116W).

NOTE:
Following my tests with real CO data, this program consumes up
to 10 times the size of the input data cube.  Especially the numpy's
std() method is not extremely memory efficient -- takes about 6 times
the size of the array.  If you provide -dTleaf and -Tcutoff on the
command line, the memory-hungry numpy routines won't be called and
the memory usage should stay below 5 times the size of your input
data cube.
"""

import sys
import os
import argparse
import numpy as np
import pyfits

__version__ = "0.1-devel"

# ============
# Main program
# ============

def main(argv=None):
    
    try:
        
        # parse arguments: if argv is None, arguments from sys.argv will be parsed
        options = parse_args(argv)
        
        # load FITS datacube
        ifits_header, idata = load_ifits(options.ifits)
        
        # derive options not set by args parser
        options = none_to_defaults(options, idata)
        
        # Init clumps mask: pixels labeled with corresponding clump numbers
        clmask = np.empty(idata.shape, dtype="int32")
        clmask[:] = -1
        """
        NOTE: dtype will be reviewed -- and changed to a smaller int sufficient for
              the number of clumps found -- before writing output FITS.
        NOTE: Clumps are first numbered from 0 on, -1 meaning no clump owning the pixel.
              Before writing the output FITS and TXT file, during renumbering, clumps will
              be numbered starting from 1, 0 meaning no clump owning the pixel.  This
              second (final) numbering will be stored in clumps' attribute final_ncl.
        """
        
        # init list of clumps
        clumps = []
        
        # Find all clumps
        print "Finding clumps."
        find_all_clumps(idata, clmask, clumps, options)
        
        # merge small clumps
        print "Merging small clumps."
        merge_small_clumps(clumps, options.Npxmin)
        
        # renumber clumps and clmask
        print "Renumbering clumps."
        renumber_clumps(clumps, options.Npxmin)
        final_clumps_count = renumber_clumps.last_new_ncl
        renumber_clmask(clmask, clumps)
        print "{0} clumps found.".format(final_clumps_count)
        """
        NOTE: Clumps have now set their final_ncl attribute, which goes from 1 on.
        NOTE: Clumps with Npx < Npxmin have their final_ncl set to 0.
        """
        
        # write clmask to output FITS
        print "Writing output FITS."
        write_ofits(options.ofits, clmask, final_clumps_count)
        
        # write clumps to output text file
        print "Writing output text file."
        write_otext(options.otext, clumps, options)
        
    except (IOError, Error) as err:
        print >>sys.stderr, err
        return 1


# =========
# Functions
# =========

def parse_args(argv=None):
    """Parses arguments using argparse module."""
    
    # setup parser
    parser = argparse.ArgumentParser(description="Identify clumps within a 3D FITS datacube.")
    parser.add_argument("-version", action="version", version=__version__)
    parser.add_argument("-ifits", required=True, help="FITS file where to search for clumps.")
    parser.add_argument("-dTleaf", type=float, help="Minimal depth of a valley separating adjacent clumps. (default: 3*sig_noise)")
    parser.add_argument("-Tcutoff", type=float, help="Minimal data value to consider. Pixels with lower values won't be processed. (default: 3*sig_noise)")
    parser.add_argument("-Npxmin", type=int, default=5, help="Minimal size of clumps in pixels. (default: %(default)s)")
    parser.add_argument("-ofits", help="FITS file where the found clumps will be saved. If exists, will be overwritten. (default: IFITS with modified extension '.clumps.fits')")
    parser.add_argument("-otext", help="Text file where the found clumps will be saved in human-readable form. If exists, will be overwritten. (default: IFITS with modified extension '.clumps.txt')")
    
    # parse args
    args = parser.parse_args(args=argv)
    
    # return namespace with parsed arguments
    return args


def load_ifits(ifits):
    """Loads and preprocess input FITS data.
    
    Returns ifits header and preprocessed idata.
    """
    
    # Load first HDU (HDU = header data unit)
    hdulist = pyfits.open(ifits)
    idata   = hdulist[0].data
    iheader = hdulist[0].header
    
    # Check if idata is 3D, i.e. has exactly 3 dimenstions.
    if idata.ndim != 3:
        raise Error("Input FITS must contain 3D data (in the first HDU), found {0}-dimensional data.".format(idata.ndim))
    
    # Add boundary to idata.
    idata2 = np.empty(np.array(idata.shape)+2, dtype=idata.dtype)
    idata2[:] = -np.inf
    idata2[1:-1, 1:-1, 1:-1] = idata
    idata = idata2
    """
    The boundary is a margin around the original data cube with values
    set to -inf.  The boundary pixels will be sorted last and the loop over
    them is expected to terminate before reaching them.  This ensures that
    all pixels with values > -inf will have their neighbours defined without
    the fear of IndexError.
    """
    
    return iheader, idata


def none_to_defaults(options, idata):
    """Derives defaults for None options.
    
    This function derives default values for options which are not set, presumably
    because they were not set by the user at the command line, and could not be set
    as simple defaults by the argument parser.  These options/arguments include:
     
     - ofits   (derived from ifits)
     - otext   (derived from ifits)
     - dTleaf  (derived from idata)
     - Tcutoff (derived from idata)
    
    Returns: updated options namespace
    """
    
    assert hasattr(options, "ifits")
    assert hasattr(options, "ofits")
    assert hasattr(options, "otext")
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    assert idata.ndim == 3
    
    new_options = options
    
    # ofits --> ifits with modified extension ".clumps.fits"
    if new_options.ofits is None:
        if   new_options.ifits[-5:] == ".fits": new_options.ofits = new_options.ifits[:-5]+".clumps.fits"
        elif new_options.ifits[-4:] == ".fit":  new_options.ofits = new_options.ifits[:-4]+".clumps.fit"
        else:                                   new_options.ofits = new_options.ifits+".clumps.fits"

    # otext --> ifits with modified extension ".clumps.txt"
    if new_options.otext is None:
        if   new_options.ifits[-5:] == ".fits": new_options.otext = new_options.ifits[:-5]+".clumps.txt"
        elif new_options.ifits[-4:] == ".fit":  new_options.otext = new_options.ifits[:-4]+".clumps.txt"
        else:                                   new_options.otext = new_options.ifits+".clumps.txt"

    # dTleaf and/or Tcutoff --> 3 * sig_noise
    if (new_options.dTleaf is None) or (new_options.Tcutoff is None):
        
        print "dTleaf and/or Tcutoff not set. Estimating from input data."
        
        # compute data mean and std
        valid = idata.view(np.ma.MaskedArray)
        valid.mask = ~np.isfinite(idata)
        mean_data = valid.mean()
        std_data = valid.std() # WARNING: Takes memory of ~ 6 times idata size.
        del valid
        
        # compute noise mean and std
        noise = idata.view(np.ma.MaskedArray)
        noise.mask = (~np.isfinite(idata)) | (idata > 3.*std_data)
        mean_noise = noise.mean()
        std_noise = noise.std() # WARNING: Takes memory of ~ 6 times idata size.
        del noise
        
        # set dTleaf
        if new_options.dTleaf is None:
            print "Setting dTleaf to {dTleaf} (= 3*std_noise = 3*{std_noise})".format(dTleaf=3.*std_noise, std_noise=std_noise)
            new_options.dTleaf = 3.*std_noise
        
        # set Tcutoff
        if new_options.Tcutoff is None:
            print "Setting Tcutoff to {Tcutoff} (= 3*std_noise = 3*{std_noise})".format(Tcutoff=3.*std_noise, std_noise=std_noise)
            new_options.Tcutoff = 3.*std_noise
    
    return new_options


def find_all_clumps(idata, clmask, clumps, options):
    """
    Finds all clumps in data cube.
    
    All means that this function will find also small clumps -- as small
    as one pixel.  These small clumps are expected to be merged or
    deleted later by other routines.
    
    Arguments:
    idata -- input 3D data cube
    clmask -- clump mask (3D array of integers)
    clumps -- list of found clumps
    options -- namespace with additional options
    """
    
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    assert idata.ndim == 3
    assert clmask.ndim == 3
    assert clmask.shape == idata.shape
    
    # Sort keys of idata array (1-D flattened keys)
    skeys1 = idata.argsort(axis=None)[::-1]

    # Find clumps -- loop over sorted keys of pixels starting at maximum
    ncl = -1 # current clump label/index
    assert options.Tcutoff > -np.inf, "Tcutoff must be > -inf." # -inf shouldn't parse through CLI
    for key1 in skeys1:
        
        # derive key3 (3-D key)
        key3 = np.unravel_index(key1, idata.shape)
        
        # get data value
        dval = idata[key3]
        
        # continue if NAN
        if dval is np.nan: continue
        
        # terminate if dval < Tcutoff
        if dval < options.Tcutoff: break
        
        # initialize pixel
        px = Pixel(key3, idata, clmask, clumps)
        
        # find neighbours (clumps touching at this pixel)
        neighbours = px.get_neighbours()
        
        # no neighbour --> new clump
        if len(neighbours) == 0:
            ncl += 1
            clumps += [Clump(ncl, px)]
            clmask[key3] = ncl
        
        # one neighbour --> add pixel to it
        elif len(neighbours) == 1:
            clmask[key3] = neighbours[0].ncl
            px.addto(neighbours[0])
        
        # more clumps --> merge/connect them
        else: # len(neighbours) > 1
            """
            There are two things to do, now:
             
              (1) Add the pixel to a clump.
              (2) Update properties of the neighbouring clumps:
                    (a) Merge too short clumps.
                    (b) Update touching lists.
                    (c) Connect grandparents.
            """
            
            # (1) Add the pixel to a clump.
            """
            Add the pixel to the nearest clump which will not be merged, i.e. to the nearest
            clump which has leaf > dTleaf, or to the tallest clump, if no clump with high
            enough leaf exists.
            """
            # find the merger (clump to which the pixel will be added)
            merger = neighbours[0] # start with tallest neighbour (neighbours are sorted)
            dist2_min = px.dist2(merger)
            for neighbour in neighbours[1:]:
                if neighbour.dpeak - px.dval < options.dTleaf:
                    break # stop search as soon as a too short clump is hit
                else:
                    dist2 = px.dist2(neighbour)
                    if dist2 < dist2_min:
                        dist2_min = dist2
                        merger = neighbour
            
            # add pixel to merger
            clmask[key3] = merger.ncl
            px.addto(merger)
            
            # (2a) Merge too short clumps.
            for neighbour in neighbours[1:]:
                if neighbour.dpeak - px.dval < options.dTleaf:
                    neighbour.parent = merger
                    neighbour.merge_to_parent()
                    """
                    NOTE: If neighbour had a parent, it wouldn't pass the IF above, since it
                    would have been already merged.
                    """
            
            # (2b) Update touching lists.
            for i, neighbour in enumerate(neighbours):
                for other_neighbour in neighbours[:i]:
                    neighbour.update_touching(other_neighbour, px.dval)
                    other_neighbour.update_touching(neighbour, px.dval)
            
            # (2c) Connect grandparents.
            gps = px.get_grandparents(neighbours) # (grandparents are sorted)
            for i in xrange(1, len(gps)):
                gp = gps[i]
                gp.parent = gps[0]
                gp.dist2_min = gp.dist2(gp.parent)
                for j in xrange(i):
                    dist2 = gp.dist2(gps[j])
                    if dist2 < gp.dist2_min:
                        gp.parent = gps[j]
                        gp.dist2_min = dist2


def merge_small_clumps(clumps, Npxmin):
    """Merges clumps with too little pixels.
    
    Clumps with too little pixels (< Npxmin) will be merged to their parents.
    Clumps will also be merged, if their parents have too little pixels.
    Merging starts from "bottom", i.e. clumps with the lowest dpeak values will
    be processed first.
    """
    
    for clump in (clumps[i] for i in xrange(len(clumps)-1, -1, -1)):
        # already merged --> skip
        if clump.merges:
            continue
        
        # solitary/orphan clump --> skip
        elif clump.parent is clump:
            continue
        
        # too small clump --> merge to its parent
        elif clump.Npx < Npxmin:
            clump.merge_to_parent()
        
        # too small parent --> merge clump to it
        elif clump.parent.get_merger().Npx < Npxmin:
            clump.merge_to_parent()


def renumber_clumps(clumps, Npxmin):
    """Renumbers clumps taking into account mergers and Npxmin limit.
    
    Sets clumps' final_ncl so that:
      (1) The numbering starts from 1.
      (2) Merged clumps are renumbered according to the final_ncl of the clump to which they merge.
          This is consistent, since clumps are expected to merge only to clumps with lower ncl, which
          are renumbered prior to the merging clump.
      (3) Solitary clumps with too little pixels (< Npxmin) are "deleted":  final_ncl is set to 0.
    
    The final count of clumps after renumbering can be retrieved from outside of this function
    through the function's attribute 'last_new_ncl'.
    """
    
    new_ncl = 0
    for clump in clumps:
        # clump merges --> use final_ncl of get_merger() clump
        if clump.merges:
            exp_clump = clump.get_merger()
            assert exp_clump.ncl < clump.ncl, "Clumps should merge only to clumps with lower ncl."
            clump.final_ncl = exp_clump.final_ncl
        
        # too small clump --> "delete"/renumber to 0
        elif clump.Npx < Npxmin:
            clump.final_ncl = 0
        
        # assign new clump number
        else:
            new_ncl += 1
            clump.final_ncl = new_ncl
    
    # save the last new_ncl (for retrieving from outside)
    renumber_clumps.last_new_ncl = new_ncl


def renumber_clmask(clmask, clumps):
    """Renumbers clmask according to clumps' final_ncl."""
    
    for ijk, ncl in np.ndenumerate(clmask):
        try:
            clmask[ijk] = clumps[ncl].final_ncl
        except KeyError:
            clmask[ijk] = 0


def write_ofits(ofits, clmask, final_clumps_count):
    """Writes clmask to output FITS file.
    
    If ofits exists it will be overwritten.
    """

    # reduce size of clmask if possible
    if final_clumps_count <= np.iinfo("uint8").max:
        clmask = clmask.astype("uint8")
    elif final_clumps_count <= np.iinfo("int16").max:
        clmask = clmask.astype("int16")

    # write FITS
    if os.path.exists(ofits): os.remove(ofits)
    pyfits.writeto(ofits, clmask[1:-1,1:-1,1:-1]) # strip clmask's boundary


def write_otext(otext, clumps, options):
    """Write clumps to output text file.
    
    The output text file is an ASCII file compatible with the original dendrofind
    textual output and as such can be used by supportive dendrofind scripts, e.g.
    df_dendrogram.py for plotting the dendrogram.
    
    If otext file exists, it will be overwritten.
    """
    
    assert hasattr(options, "Tcutoff")
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Npxmin")
    
    with open(otext, "w") as f:
        f.write("# Nlevels = 1000 Tcutoff = {options.Tcutoff} dTleaf = {options.dTleaf} Npxmin = {options.Npxmin}\n".format(options=options))
        """
        The Nlevels value is set to 1000 only for the output to be compatible with the original dendrofind's
        textual output.  It has no real meaning for df2, since df2 doesn't use Nlevels parameter.
        """

        # output clumps
        for clump in clumps:
            # print only clumps which were not deleted (final_ncl > 0) or merged
            if (not clump.merges) and (clump.final_ncl > 0):
                f.write(str(clump))
                f.write("\n")


# =======
# Classes
# =======

class Error(Exception):
    """Standard non-specific runtime error.
    
    This exception is intended as a correct and single exit point for main()
    in the case an anticipated error occurs.
    """
    
    def __init__(self, msg):
        self.msg = msg
    
    def __str__(self):
        return "Error: {0}".format(self.msg)


class Clump(object):
    """Clump found within the data cube."""
    
    def __init__(self, ncl, px):
        """Build a new clump.
        
        ncl -- number/label for this new clump
        px  -- the first pixel of the clump (peak)
        """
        
        self.ncl = ncl         # label
        self.final_ncl = None  # final clump label to be set during renumbering
        self.Npx = 1           # number of pixels (incl. pixels of clumps which merge to this clump)
        self.pixels = [[px.ijk, px.dval]] # list of pixels belonging to the clump: [[(x,y,z), dval], ...]
        self.parent = self     # parent to which the clump is connected (the nearest clump with
                               # higher dpeak touching this clump)
        self.merges = False    # whether the clump merges to its parent
        self.touching = {}     # dictionary of other clumps which touch this one: {clump_reference: dval_at_which_they_touch}
        self.xyz = px.xyz      # (x,y,z) coordinates: the weighted average with the weight equal to pixels'
                               # data values. Note xyz changes as new pixels are being added to the clump.
        self.wxyz = px.dval    # weight of clump's xyz coordinates
        self.dpeak = px.dval   # peak data value
        self.sumd = px.dval    # sum of clumps data values (sum of all pixels' dval)
    
    def get_merger(self):
        """Returns clump to which this clump merges or self."""
        return self.parent.get_merger() if self.merges else self
    
    def get_grandparent(self):
        """Returns parent's parent's... parent or self."""
        return self.parent.get_grandparent() if self.parent is not self else self
    
    def dist2(self, other):
        """Returns square of the distance of the clump to other clump or pixel."""
        return sum( (self.xyz-other.xyz)**2 )
    
    def merge_to_parent(self):
        """Merges clump to its parent.
        
        The parent is expanded, i.e., the clump is merged to self.parent.get_merger().
        """
        
        assert not self.merges, "Attempt to merge already merged clump."
        
        # get the clump to which merge
        merger = self.parent.get_merger()
        
        assert merger is not self, "Attempt to merge clump to itself."
        
        # ------ merge ------
        # update pixels
        merger.Npx += self.Npx
        merger.pixels.extend(self.pixels)
        
        # update xyz
        merger.xyz = (self.wxyz*self.xyz + merger.wxyz*merger.xyz)/(self.wxyz + merger.wxyz)
        merger.wxyz += self.wxyz
        merger.sumd += self.sumd
        
        # update touching
        for clump, touching_at_dval in self.touching.iteritems():
            merger.update_touching(clump, touching_at_dval)
        
        # set merges tag
        self.merges = True
     
    def update_touching(self, other, dval):
        """Updates dict of clumps which touch this clump.
        
        Clump "other" is expanded to other.get_merger().
        dval is the data value of the pixel which connects clumps "self" and "other".
        """
        
        # expand other
        exp_other = other.get_merger()
        
        # update touching dict
        if exp_other is not self:
            if (exp_other not in self.touching) or (dval > self.touching[exp_other]):
                self.touching.update({exp_other: dval})
    
    def compact_touching(self):
        """Compacts touching dict.
        
        (1) Removes references to deleted clumps (with final_ncl == 0).
        (2) Ensures the touching dict is unique (solving mergers).  If more references to
            the same touched clump are found, sets the touching_at_dval to the highest value.
        """
        
        new_touching = {}
        for clump, touching_at_dval in self.touching.iteritems():
            exp_clump = clump.get_merger() # expand touched clump
            if exp_clump.final_ncl == 0: continue # ditch deleted clumps (with final_ncl == 0)
            if (exp_clump not in new_touching) or (touching_at_dval > new_touching[exp_clump]):
                new_touching.update({exp_clump: touching_at_dval})
        self.touching = new_touching
    
    def get_connected(self):
        """Returns clumps connected to this clump.
        
        Connected clumps are all the clumps which either touch a given clump
        directly or indirectly through other connected clumps.  This structure
        is used for building the dendrogram.  In other words, connected clumps make
        up a graph data structure.  We now want to find all the clumps (nodes)
        connected to a given clump (self) -- i.e. discover the whole graph.
        
        Returns: connected -- dict of connected clumps in form {clump: connects_at_dval, ...},
                              where connects_at_dval is the data value at which the clump connects
        """
        
        # Init the queue of clumps to explore, start with self.get_merger().
        # Queue format: [[clump, dval], ...], where dval is the data value at which the clump connects
        # NOTE: only expanded clumps are expected in the queue
        queue = [[self.get_merger(), self.get_merger().dpeak]]
        
        # Init the connected dict
        # NOTE: only expanded clumps are expected in the dict
        connected = dict(queue)
        
        # find all connected clumps (discover the whole graph, incl. clump self)
        while queue:
            next_in_queue = queue.pop() # LIFO queue --> depth-first traversal
            focused_clump = next_in_queue[0]
            focused_dval = next_in_queue[1]
            
            assert not focused_clump.merges, "Only expanded clumps are expected in the queue."
            
            for child_clump, child_dval in focused_clump.touching.iteritems():
                exp_child_clump = child_clump.get_merger() # expand clump
                
                # get the minimal data value along the path
                min_dval = min(focused_dval, child_dval)
                
                # newly discovered clump
                if exp_child_clump not in connected:
                    queue.append([exp_child_clump, min_dval])
                    connected.update({exp_child_clump: min_dval})
                # rediscovered clump
                else:
                    # update if found a better/"higher" path (with greater minimal dval along it)
                    if min_dval > connected[exp_child_clump]:
                        queue.append([exp_child_clump, min_dval])
                        connected[exp_child_clump] = min_dval
        
        # remove self.get_merger() from connected
        del connected[self.get_merger()]
        
        return connected
    
    def __str__(self):
        """Returns textual representation of the clump.
        
        The output text is compatible with the original dendrofind's textual output and should
        work with the script for plotting dendrograms (df_dendrogram.py).
        """
        
        # compact touching clumps dict first
        self.compact_touching()
        
        # sort touching clumps (order by dval_at_which_they_touch, ncl)
        touching = list(self.touching.iteritems())
        touching.sort(key=lambda x: (-x[1], x[0].final_ncl))
        
        # get connected clumps, then sort them (by dval_at_which_they_connect, ncl)
        connected = self.get_connected()
        connected = list(connected.iteritems())
        connected.sort(key=lambda x: (-x[1], x[0].final_ncl))
        
        # sort list of pixels (order by dval, k, j, i)
        self.pixels.sort(key=lambda x: (-x[1], x[0][2], x[0][1], x[0][0]))
        
        # generate text_repr to be returned
        text_repr  = "clump: {final_ncl}\n".format(final_ncl=self.final_ncl)
        text_repr += "  Npx: {Npx}\n".format(Npx=self.Npx)
        text_repr += "  Tmax: {dpeak}\n".format(dpeak=self.dpeak)
        text_repr += "  state: independent\n"
        text_repr += "  Ntouching: {Ntouching}\n".format(Ntouching=len(touching))
        text_repr += "  Nconnected: {Nconnected}\n".format(Nconnected=len(connected))
        text_repr += "  pixels:\n"
        for px in self.pixels:
            text_repr += "    {ijk[2]:>3d} {ijk[1]:>3d} {ijk[0]:>3d} {dval}\n".format(ijk=px[0], dval=px[1])
            """
            NOTE: PyFITS reverses the order of coordinates, therefore we output 2-1-0.
            NOTE: Coordinates in FITS start from 1 and since arrays clmask and idata in this code
                  added an one-pixel border around the original data, these two shifts cancel each
                  other out and we can output ijk directly.
            """
        text_repr += "  touching:\n"
        for cl in touching:
            text_repr += "    {final_ncl:>3d} {dval}\n".format(final_ncl=cl[0].final_ncl, dval=cl[1])
        text_repr += "  connected:\n"
        for cl in connected:
            text_repr += "    {final_ncl:>3d} {dval}\n".format(final_ncl=cl[0].final_ncl, dval=cl[1])
        
        return text_repr


class Pixel(object):
    """Pixel within the data cube."""
    
    # relative map for neighbouring pixels
    neigh_map = np.array([ (0,0,+1),
                           (0,0,-1),
                           (0,+1,0),
                           (0,-1,0),
                           (+1,0,0),
                           (-1,0,0) ], dtype=int)
    
    def __init__(self, ijk, idata, clmask, clumps):
        self.ijk = np.array(ijk, dtype=int)   # (i,j,k) coordinates
        self.xyz = np.array(ijk, dtype=float) # (x,y,z) coordinates
        self.idata = idata                    # input data cube
        self.clmask = clmask                  # pixels' mask (to which clump the pixels belong)
        self.clumps = clumps                  # list of clumps
        self.dval = idata[tuple(self.ijk)]    # data value
    
    def get_neighbours(self):
        """Finds neighbours touching at this pixel.
        
        The list of neighbours is sorted: clumps with lower ncl first.
        """
        
        neighbours = []
        for shift in self.neigh_map:
            ncl = self.clmask[tuple(self.ijk + shift)]
            if ncl > -1:
                neighbour = self.clumps[ncl].get_merger()
                if neighbour not in neighbours:
                    neighbours.append(neighbour)
        neighbours.sort(key=lambda clump: clump.ncl) # clumps with lower ncl first
        
        return neighbours
    
    def get_grandparents(self, neighbours=None):
        """Finds grandparents (parent's parent's... parent) among neighbours.
        
        The grandparents are searched for in the list of neighbours.  If the list is not
        provided, it will be taken from self.get_neighbours().
        
        The list of grandparents is sorted: clumps with lower ncl first.
        """
        
        if neighbours is None: neighbours = self.get_neighbours()
        
        grandparents = []
        for neighbour in neighbours:
            grandparent = neighbour.get_grandparent()
            if grandparent not in grandparents:
                grandparents.append(grandparent)
        grandparents.sort(key=lambda clump: clump.ncl) # clumps with lower ncl first
        
        return grandparents
    
    def dist2(self, other):
        """Returns square of the distance of the pixel to other pixel or clump."""
        return sum( (self.xyz-other.xyz)**2 )
    
    def addto(self, clump):
        """Adds pixel to clump."""
        clump.Npx += 1
        clump.pixels.append([self.ijk, self.dval])
        clump.xyz = (self.dval*self.xyz + clump.wxyz*clump.xyz)/(self.dval + clump.wxyz) # use dval as the weight
        clump.wxyz += self.dval
        clump.sumd += self.dval


# ====================
# Execute main program
# ====================

if __name__ == "__main__":
    sys.exit(main())
