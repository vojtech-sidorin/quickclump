#!/usr/bin/env python
# clumpit - identify clumps within a 3D FITS datacube
#
# Copyright 2015 Vojtech Sidorin <vojtech.sidorin@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Identify clumps within a 3D FITS datacube.

QUICK START GUIDE
=================

To find clumps in your FITS data cube, type

    $ python clumpit.py my_datacube.fits

To show usage help, type

    $ python clumpit.py -h

To run in the interactive mode, type

     $ python
     >>> import clumpit
     >>> clumpit.main(["my_datacube.fits"])

DESCRIPTION
===========

This program is an improved implementation of DENDROFIND(1) -- a clump-
finding algorithm inspired by Clumpfind(2).  DENDROFIND was originally
conceived by Richard Wunsch, who also published its first
implementation in Python, later rewritten in C.  Compared to the
original implementation, clumpit uses different data structures and
doesn't need parameter Nlevels.  clumpit is also faster (about 50 000
times) and scales linearly with the data cube volume (number of pixels).

(1) See <http://galaxy.asu.cas.cz/~richard/dendrofind/> for a
    description of the original DENDROFIND algorithm.  The first
    practical use together with another description was published by
    Wunsch et al. (2012), see
    <http://adsabs.harvard.edu/abs/2012A%26A...539A.116W>.
(2) See <http://www.ifa.hawaii.edu/users/jpw/clumpfind.shtml> or
    <http://adsabs.harvard.edu/abs/1994ApJ...428..693W>.

NOTES
=====

Following my tests with real CO data, this program consumes up to
10 times the size of the input data cube.  Numpy's std() method is
especially eager for memory and takes about 6 times the size of the
array (input data cube).  However, if you provide the parameters
--dTleaf and --Tcutoff at the command-line, the memory-hungry numpy
routines won't be called and the memory usage should stay below 5 times
the size of your input data cube.

Tested in Python 2.7 and 3.4.
"""

import sys
import os
import argparse
import datetime

import numpy as np
# Import FITS IO.
# NOTE: PyFITS was merged into Astropy.
try:
    from astropy.io import fits
except ImportError:
    try:
        import pyfits as fits
    except ImportError:
        sys.exit("Error: Cannot find any supported FITS IO package.  "
                 "Do you have installed 'astropy' or 'pyfits'?")

__version__ = "1.3-2"



def main(argv=None):
    try:
        _main(argv=argv)
    except (IOError, InputDataError, OutOfBoundsError) as e:
        sys.stderr.write("{0}: {1}\n".format(e.__class__.__name__, str(e)))
        return 1

def _main(argv=None):
    # Parse arguments: if argv is None, arguments from sys.argv will be
    # used automatically.
    options = parse_args(argv)

    # Load the input data (a FITS datacube).
    idata = load_idata(options.ifits)

    # Set options that were not set by the args parser.
    options = set_defaults(options, idata)

    check_options(options)

    # Initialise the clumps mask: pixels labeled with the number of the
    # corresponding clump.
    clmask = np.empty(idata.shape, dtype="int32")
    clmask[:] = -1
    # NOTE: dtype will be reviewed later, before saving the output into a
    # FITS file, and changed to a smaller int sufficient for holding the
    # number of the found clumps.
    # NOTE: Initially, the clumps are numbered from 0 on, -1 meaning no
    # clump owning the pixel.  The final numbering, stored in clumps'
    # atribute final_ncl, will start from 1 with 0 meaning no clumps owning
    # the pixel.

    # init list of clumps
    clumps = []

    if options.verbose > 0:
        print("Finding clumps.")
    find_all_clumps(idata, clmask, clumps, options)

    if options.verbose > 0:
        print("Merging small clumps.")
    merge_small_clumps(clumps, options.Npxmin)

    if options.verbose > 0:
        print("Renumbering clumps.")
    final_clumps_count = renumber_clumps(clumps, options.Npxmin)
    renumber_clmask(clmask, clumps)
    if options.verbose >= 0:
        print("{N} clumps found.".format(N=final_clumps_count))
    # NOTE: The clumps have now set their final labels/numbers, which are
    # stored in attribute final_ncl.
    # NOTE: Too small clumps, those with Npx < Npxmin, have set their
    # final_ncl to 0.

    if options.ofits.strip().upper() != "NONE":
        if options.verbose > 0:
            print("Writing output FITS.")
        write_ofits(options.ofits, clmask, final_clumps_count, options)

    if options.otext.strip().upper() != "NONE":
        if options.verbose > 0:
            print("Writing output text file.")
        write_otext(options.otext, clumps, options)

def parse_args(argv=None):
    """Parse arguments with argparse."""
    parser = argparse.ArgumentParser(description="Identifies clumps within a "
                                     "3D FITS datacube.")
    parser.add_argument("ifits", help="FITS file where to search for clumps.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--dTleaf", type=float, help="Minimal depth of a "
                        "valley separating adjacent clumps.  Must be > 0.  "
                        "(default: 3*sig_noise)")
    parser.add_argument("--Tcutoff", type=float, help="Minimal data value to "
                        "consider.  Pixels with lower values won't be "
                        "processed.  Must be > 0.  (default: 3*sig_noise)")
    parser.add_argument("--Npxmin", type=int, default=5, help="Minimal size "
                        "of a clump in pixels.  (default: %(default)s)")
    parser.add_argument("--ofits", help="FITS file where the found clumps "
                        "will be saved.  If OFITS exists, it will be "
                        "overwritten.  If set to 'None' (case doesn't "
                        "matter), OFITS file won't be written.  "
                        "(default: IFITS with modified extension "
                        "'.clumps.fits')")
    parser.add_argument("--otext", help="Text file where the found clumps "
                        "will be saved in a human-readable form.  If OTEXT "
                        "exists, it will be overwritten.  If set to 'None' "
                        "(case insensitive), OTEXT file won't be written.  "
                        "This will speed up the program's execution.  On the "
                        "other hand, the OTEXT file is needed for the "
                        "construction of a dendrogram.  "
                        "(default: IFITS with modified extension "
                        "'.clumps.txt')")
    parser.add_argument("--verbose", "-v", action="count", default=0,
                        help="Increase verbosity.")
    parser.add_argument("--silent", dest="verbose", action="store_const",
                        const=-1, help="Suppress output to stdout.  (Set "
                        "verbosity to a minimum.)")
    args = parser.parse_args(argv)
    return args

def load_idata(ifits):

    """Load and preprocess input FITS data."""

    # Load the first HDU from the FITS (HDU = header data unit).
    with fits.open(ifits) as f:
        idata = f[0].data  # f[0] == the first HDU in the file.

    # Check if idata is 3D, i.e. has exactly 3 dimensions.
    if idata.ndim != 3:
        raise InputDataError("The input FITS file must contain 3D data (in "
                             "the first HDU), found {0}-dimensional data."
                             .format(idata.ndim))

    # Add boundary to idata.
    idata2 = np.empty(np.array(idata.shape)+2, dtype=idata.dtype)
    idata2[:] = -np.inf
    idata2[1:-1, 1:-1, 1:-1] = idata
    idata = idata2
    # NOTE: The boundary is a margin around the original data cube with values
    # set to -inf.  The boundary pixels will be sorted last and the loop over
    # them is expected to terminate before reaching them.  This ensures that
    # all pixels with values > -inf will have their neighbours defined without
    # the fear of IndexError.

    return idata

def set_defaults(options, idata):

    """Set default values for options.

    This function derives default values for options which are not set,
    presumably because they were not set by the user at the command line,
    and could not be set as simple defaults by the argument parser.
    These options include:

     - ofits   (derived from ifits)
     - otext   (derived from ifits)
     - dTleaf  (derived from idata)
     - Tcutoff (derived from idata)

    Return updated options namespace.
    """

    assert hasattr(options, "ifits")
    assert hasattr(options, "ofits")
    assert hasattr(options, "otext")
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    assert idata.ndim == 3

    new_options = options

    # ofits -- ifits with modified extension ".clumps.fits"
    if new_options.ofits is None:
        if new_options.ifits.endswith(".fits"):
            new_options.ofits = new_options.ifits[:-5] + ".clumps.fits"
        elif new_options.ifits.endswith(".fit"):
            new_options.ofits = new_options.ifits[:-4] + ".clumps.fit"
        else:
            new_options.ofits = new_options.ifits + ".clumps.fits"

    # otext -- ifits with modified extension ".clumps.txt"
    if new_options.otext is None:
        if new_options.ifits.endswith(".fits"):
            new_options.otext = new_options.ifits[:-5] + ".clumps.txt"
        elif new_options.ifits.endswith(".fit"):
            new_options.otext = new_options.ifits[:-4] + ".clumps.txt"
        else:
            new_options.otext = new_options.ifits + ".clumps.txt"

    # dTleaf/Tcutoff -- 3*sig_noise
    if (new_options.dTleaf is None) or (new_options.Tcutoff is None):
        if options.verbose > 0:
            print("dTleaf and/or Tcutoff not set.  Estimating from the input "
                  "data (IFITS).")

        # Compute data mean and std.
        valid = idata.view(np.ma.MaskedArray)
        valid.mask = ~np.isfinite(idata)
        # NOTE: Numpy's std() takes memory of ~6 times idata size.
        std_data = valid.std()
        del valid  # No longer needed: delete the reference

        # Compute noise mean and std.
        noise = idata.view(np.ma.MaskedArray)
        noise.mask = (~np.isfinite(idata)) | (idata > 3.*std_data)
        # NOTE: Numpy's std() takes memory of ~6 times idata size.
        std_noise = noise.std()
        del noise  # No longer needed: delete the reference

        # Check if estimation of std_noise from input data succeeded.
        if (not np.isfinite(std_noise)) or (std_noise <= 0.):
            raise OutOfBoundsError(
                "Estimation of std_noise from input data failed.  "
                "Got value '{0}'.  "
                "Is the input FITS data valid/reasonable?".format(std_noise))

        # Set dTleaf.
        if new_options.dTleaf is None:
            if options.verbose > 0:
                print("Setting dTleaf to {dTleaf} (= 3*std_noise = "
                      "3*{std_noise})"
                      .format(dTleaf=3.*std_noise, std_noise=std_noise))
            new_options.dTleaf = 3.*std_noise

        # Set Tcutoff.
        if new_options.Tcutoff is None:
            if options.verbose > 0:
                print("Setting Tcutoff to {Tcutoff} (= 3*std_noise = "
                      "3*{std_noise})"
                      .format(Tcutoff=3.*std_noise, std_noise=std_noise))
            new_options.Tcutoff = 3.*std_noise

    return new_options

def check_options(options):
    """Check values of dTleaf and Tcutoff."""
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    if not (options.dTleaf > 0.):
        raise OutOfBoundsError("'dTleaf' must be > 0.")
    if not (options.Tcutoff > 0.):
        raise OutOfBoundsError("'Tcutoff' must be > 0.")

def find_all_clumps(idata, clmask, clumps, options):

    """Find all clumps in data cube idata.

    'All' means that this function will also find small clumps --
    as small as one pixel.  The clumps that are smaller than Npxmin
    are expected to be merged or deleted later by other routines
    following this function.

    Positional arguments:

     idata   -- input 3D data cube
     clmask  -- clump mask (3D array of integers)
     clumps  -- list of found clumps
     options -- namespace with additional options

    This function updates clmask and clumps in-place.
    """

    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    assert idata.ndim == 3
    assert clmask.ndim == 3
    assert clmask.shape == idata.shape

    # Sort keys of idata array (1-D flattened keys).
    skeys1 = idata.argsort(axis=None)[::-1]

    # Find clumps -- loop over sorted keys of pixels starting at maximum.
    ncl = -1  # current clump label/index
    assert options.Tcutoff > 0.
    for key1 in skeys1:

        # Derive key3 (3-D key)
        key3 = np.unravel_index(key1, idata.shape)

        # Get data value
        dval = idata[key3]

        # Skip NANs
        if dval is np.nan:
            continue

        # Terminate if dval < Tcutoff.  Since the keys are sorted, we can
        # terminate the loop.
        if dval < options.Tcutoff:
            break

        # Initialize pixel
        px = Pixel(key3, idata, clmask, clumps)

        # Find neighbours (clumps touching at this pixel)
        neighbours = px.get_neighbours()

        if not neighbours:
            # No neighbour --> Make a new clump
            ncl += 1
            clumps.append(Clump(ncl, px))
            clmask[key3] = ncl
        elif len(neighbours) == 1:
            # One neighbour --> Add pixel to it
            clmask[key3] = neighbours[0].ncl
            px.addto(neighbours[0])
        else:
            # More neighbours --> Merge/connect them
            # NOTE: There are two things to do now:
            #
            #  (1) Add the pixel to a clump.
            #  (2) Update the properties of the neighbouring clumps:
            #       (a) Merge too short clumps.
            #       (b) Update touching lists.
            #       (c) Connect grandparents.

            # (1) Add the pixel to a clump.
            # NOTE: Add the pixel to the nearest clump which will not be
            # merged, i.e. to the nearest clump which has leaf > dTleaf, or to
            # the tallest clump if no clump with a leaf that is high enough
            # exists.

            # Find the merger (clump to which the pixel will be added)
            # Start with the tallest neighbour (neighbours are sorted)
            merger = neighbours[0]
            dist2_min = px.dist2(merger)
            for neighbour in neighbours[1:]:
                if neighbour.dpeak - px.dval < options.dTleaf:
                    break  # Stop search as soon as a too short clump is hit
                else:
                    dist2 = px.dist2(neighbour)
                    if dist2 < dist2_min:
                        dist2_min = dist2
                        merger = neighbour

            # Add pixel to merger
            clmask[key3] = merger.ncl
            px.addto(merger)

            # (2) Update the properties of the neighbouring clumps.
            # (2a) Merge too short clumps.
            for neighbour in neighbours[1:]:
                if neighbour.dpeak - px.dval < options.dTleaf:
                    neighbour.parent = merger
                    neighbour.merge_to_parent()
                    # NOTE: If the neighbour had a parent, it wouldn't pass the
                    # IF above, since it would have been already merged.

            # (2b) Update touching lists.
            for i, neighbour in enumerate(neighbours):
                for other_neighbour in neighbours[:i]:
                    neighbour.update_touching(other_neighbour, px.dval)
                    other_neighbour.update_touching(neighbour, px.dval)

            # (2c) Connect grandparents.
            # NOTE: Grandparents are sorted.
            gps = px.get_grandparents(neighbours)
            for i in range(1, len(gps)):
                gp = gps[i]
                gp.parent = gps[0]
                gp.dist2_min = gp.dist2(gp.parent)
                for j in range(i):
                    dist2 = gp.dist2(gps[j])
                    if dist2 < gp.dist2_min:
                        gp.parent = gps[j]
                        gp.dist2_min = dist2

def merge_small_clumps(clumps, Npxmin):
    """Merge clumps with too little pixels.

    Clumps with too little pixels (< Npxmin) will be merged to their
    parents.  Clumps will also be merged, if their parents have too little
    pixels.  Merging starts from the "bottom", i.e. clumps with the lowest
    dpeak values will be processed first.
    """
    for clump in reversed(clumps):
        if clump.merges:
            # Already merged --> skip
            continue
        elif clump.parent is clump:
            # Solitary/orphan clump --> skip
            continue
        elif clump.Npx < Npxmin:
            # Too small clump --> merge to its parent
            clump.merge_to_parent()
        elif clump.parent.get_merger().Npx < Npxmin:
            # Too small parent --> merge clump to it
            clump.merge_to_parent()

def renumber_clumps(clumps, Npxmin):
    """Renumber clumps taking into account mergers and Npxmin limit.

    Set clumps' final_ncl so that:

     - The numbering starts from 1.
     - Merged clumps are renumbered according to the final_ncl of the clump
       to which they merge.  This is consistent, since clumps are expected to
       merge only to clumps with a lower ncl, which are processed/renumbered
       prior to the merging clump.
     - Solitary clumps with too little pixels (< Npxmin) are "deleted",
       i.e. final_ncl is set to 0.

    Return the final count of clumps: last new_ncl.
    """
    new_ncl = 0
    for clump in clumps:
        if clump.merges:
            # Clump merges --> use final_ncl of get_merger() clump
            exp_clump = clump.get_merger()
            assert exp_clump.ncl < clump.ncl, \
                    "Clumps should merge only to clumps with lower ncl."
            clump.final_ncl = exp_clump.final_ncl
        elif clump.Npx < Npxmin:
            # Too small clump --> renumber to 0 (delete)
            clump.final_ncl = 0
        else:
            # Assign new clump number
            new_ncl += 1
            clump.final_ncl = new_ncl
    return new_ncl

def renumber_clmask(clmask, clumps):
    """Renumber clmask according to clumps' final_ncl."""
    if not clumps:
        clmask[:] = 0
    else:
        for ijk, ncl in np.ndenumerate(clmask):
            if clmask[ijk] < 0:
                clmask[ijk] = 0
            else:
                clmask[ijk] = clumps[ncl].final_ncl

def write_ofits(ofits, clmask, final_clumps_count, options):

    """Write clmask to the output FITS file (ofits).

    If the output FITS (ofits) exists it will be overwritten.
    """

    assert hasattr(options, "ifits")
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    assert hasattr(options, "Npxmin")

    # Reduce the size of clmask if possible.
    if final_clumps_count <= np.iinfo("uint8").max:
        clmask = clmask.astype("uint8")
    elif final_clumps_count <= np.iinfo("int16").max:
        clmask = clmask.astype("int16")

    # Create a new FITS HDU.  Compensate for the border.
    ohdu = fits.PrimaryHDU(clmask[1:-1, 1:-1, 1:-1])

    # Set the header.
    ohdu.header["BUNIT"] = ("Ncl", "clump number")
    ohdu.header["DATE"] = (
        datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S UTC"),
        "file creation date"
        )
    ohdu.header["COMMENT"] = ("File created by clumpit.py (v{version})."
                              .format(version=__version__))
    ohdu.header["COMMENT"] = ("Original data file: '{ifits}'"
                              .format(ifits=os.path.basename(options.ifits)))
    ohdu.header["COMMENT"] = ("clumpit was run with following parameters:")
    ohdu.header["COMMENT"] = ("  dTleaf={dTleaf}".format(dTleaf=options.dTleaf))
    ohdu.header["COMMENT"] = ("  Tcutoff={Tcutoff}"
                              .format(Tcutoff=options.Tcutoff))
    ohdu.header["COMMENT"] = ("  Npxmin={Npxmin}".format(Npxmin=options.Npxmin))
    ohdu.header["COMMENT"] = ("Total clumps found: {fcc}"
                              .format(fcc=final_clumps_count))
    ohdu.header["COMMENT"] = (
        "This FITS file is a mask for the original data file '{ifits}'. "
        "Each pixel contains an integer that corresponds to the label of the "
        "clump that owns the pixel. Pixels marked with zeroes belong to no "
        "clump.".format(ifits=os.path.basename(options.ifits))
        )

    # Write the output FITS.
    if os.path.exists(ofits):
        os.remove(ofits)
    ohdu.writeto(ofits)

def write_otext(otext, clumps, options):

    """Write clumps to the output text file.

    The output text file is compatible with DENDROFIND's textual output
    and as such can be used as an input for supportive DENDROFIND's
    scripts, e.g. df_dendrogram.py for plotting the dendrogram.

    If the output text file (otext) exists it will be overwritten.
    """

    assert hasattr(options, "Tcutoff")
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Npxmin")

    with open(otext, "w") as f:
        # Output the header line.
        # NOTE: The Nlevels value is set to 1000 only for the output to be
        # compatible with DENDROFIND.  It has no meaning for clumpit.
        f.write("# Nlevels = 1000 Tcutoff = {options.Tcutoff:.12g} dTleaf = "
                "{options.dTleaf:.12g} Npxmin = {options.Npxmin}\n"
                .format(options=options))
        # Output the clumps.
        for clump in clumps:
            # Output only clumps which were not deleted (final_ncl > 0) or
            # merged.
            if (not clump.merges) and (clump.final_ncl > 0):
                f.write(str(clump))
                f.write("\n")


class InputDataError(Exception):
    """Error in the input FITS data."""


class OutOfBoundsError(Exception):
    """Error signaling that a variable is out of expected bounds."""


class Clump(object):

    """Clump found within the data cube."""

    def __init__(self, ncl, px):
        """Build a new clump with first pixel px.

        Positional arguments:

         ncl -- number/label for this clump
         px  -- the first pixel of the clump (peak)
        """
        # Clump label
        self.ncl = ncl
        # Final clump label to be set during last renumbering.
        self.final_ncl = None
        # Number of pixels (incl. pixels of clumps which merge to this clump).
        self.Npx = 1
        # List of pixels belonging to the clump: [[(x,y,z), dval], ...]
        self.pixels = [[px.ijk, px.dval]]
        # Parent to which the clump is connected (the nearest clump with a
        # higher dpeak touching this clump).
        self.parent = self
        # Whether the clump merges to its parent.
        self.merges = False
        # Other clumps which touch this one.
        # {clump_reference: dval_at_which_they_touch}
        self.touching = {}
        # (x,y,z) coordinates of the clump:  The weighted average with the
        # weight equal to the clump's pixels data values.  Note the xyz changes
        # as new pixels are being added to the clump.
        self.xyz = px.xyz
        # The weight of the clump's xyz coordinates.
        self.wxyz = px.dval
        # The peak data value.
        self.dpeak = px.dval
        # The sum of the clump's data values (sum of all pixels' dval).
        self.sumd = px.dval

    def get_merger(self):
        """Return clump to which this clump merges or self."""
        return self.parent.get_merger() if self.merges else self

    def get_grandparent(self):
        """Return parent's parent's... parent or self."""
        if self.parent is not self:
            return self.parent.get_grandparent()
        else:
            return self

    def dist2(self, other):
        """Return square of the distance to the other clump or pixel."""
        return ((self.xyz-other.xyz)**2).sum()

    def merge_to_parent(self):

        """Merge clump to its parent.

        The parent is expanded, i.e. the clump is merged to
        self.parent.get_merger().
        """

        assert not self.merges, "Attempt to merge already merged clump."

        # Get the clump to which merge.
        merger = self.parent.get_merger()

        assert merger is not self, "Attempt to merge clump to itself."

        # Update pixels
        merger.Npx += self.Npx
        merger.pixels.extend(self.pixels)

        # Update xyz
        merger.xyz = ((self.wxyz*self.xyz + merger.wxyz*merger.xyz)/
                      (self.wxyz + merger.wxyz))
        merger.wxyz += self.wxyz
        merger.sumd += self.sumd

        # Update touching
        for clump, touching_at_dval in self.touching.items():
            merger.update_touching(clump, touching_at_dval)

        # Set merges tag
        self.merges = True

    def update_touching(self, other, dval):
        """Add clump other to the dict of touching clumps.

        Positional arguments:

         other -- The clump which touches this one.  Note that clump other is
                  expanded to other.get_merger().
         dval  -- The data value of the pixel which connects the two clumps.
        """
        exp_other = other.get_merger()
        if exp_other is not self:
            if ((exp_other not in self.touching) or
                    (dval > self.touching[exp_other])):
                self.touching.update({exp_other: dval})

    def compact_touching(self):
        """Compact the touching dict.

        (1) Remove references to deleted clumps (with final_ncl == 0).
        (2) Ensure the touching dict is unique (solving mergers).  If more
            references to the same clump are found, sets the touching_at_dval
            to the highest value.
        """
        new_touching = {}
        for clump, touching_at_dval in self.touching.items():
             # Expand the touched clump.
            exp_clump = clump.get_merger()
            if exp_clump.final_ncl == 0:
                # Ditch deleted clumps (with final_ncl == 0).
                continue
            if ((exp_clump not in new_touching) or
                    (touching_at_dval > new_touching[exp_clump])):
                new_touching.update({exp_clump: touching_at_dval})
        self.touching = new_touching

    def get_connected(self):

        """Return clumps connected to this clump.

        Connected clumps are all the clumps which either touch this clump
        directly or indirectly through other connected clumps.  This structure
        is used for building the dendrogram.  In other words, connected clumps
        make up a graph data structure.  We now want to find all the clumps
        (nodes) connected to clump self -- i.e. to discover the whole graph.

        Return connected -- dict of connected clumps in the form of
        {clump: connects_at_dval, ...}, where connects_at_dval is the data
        value at which the clump connects.
        """

        # Init the queue of clumps to explore; start with self.get_merger().
        # Queue format: [[clump, dval], ...], where dval is the data value at
        # which the clump connects.
        # NOTE: Only expanded clumps are expected in the queue.
        queue = [[self.get_merger(), self.get_merger().dpeak]]

        # Init the connected dict.
        # NOTE: Only expanded clumps are expected in the dict
        connected = dict(queue)

        # Find all connected clumps (discover the whole graph, incl. clump
        # self).
        while queue:
            # LIFO queue --> depth-first traversal
            next_in_queue = queue.pop()
            focused_clump = next_in_queue[0]
            focused_dval = next_in_queue[1]
            assert not focused_clump.merges, \
                "Only expanded clumps are expected in the queue."
            for child_clump, child_dval in focused_clump.touching.items():
                # Expand the clump.
                exp_child_clump = child_clump.get_merger()
                # Get the minimal data value along the path.
                min_dval = min(focused_dval, child_dval)
                if exp_child_clump not in connected:
                    # Newly discovered clump
                    queue.append([exp_child_clump, min_dval])
                    connected.update({exp_child_clump: min_dval})
                else:
                    # Rediscovered clump; update if found a better/"higher"
                    # path (with greater minimal dval along it).
                    if min_dval > connected[exp_child_clump]:
                        queue.append([exp_child_clump, min_dval])
                        connected[exp_child_clump] = min_dval

        # Remove self.get_merger() from connected.
        del connected[self.get_merger()]

        return connected

    def __str__(self):

        """Return textual representation of the clump.

        The output text is compatible with the original dendrofind's textual
        output and should work with the script for plotting dendrograms
        (df_dendrogram.py).
        """

        # Compact touching clumps dict first.
        self.compact_touching()

        # Sort touching clumps (order by dval_at_which_they_touch, ncl).
        touching = list(self.touching.items())
        touching.sort(key=lambda x: (-x[1], x[0].final_ncl))

        # Get connected clumps.
        # Then sort them (by dval_at_which_they_connect, ncl).
        connected = self.get_connected()
        connected = list(connected.items())
        connected.sort(key=lambda x: (-x[1], x[0].final_ncl))

        # Sort list of pixels (order by dval, k, j, i).
        self.pixels.sort(key=lambda x: (-x[1], x[0][2], x[0][1], x[0][0]))

        # Generate str_ to be returned.
        str_ = ["clump: {final_ncl}\n"
                "  Npx: {Npx}\n"
                "  Tmax: {Tmax:.12g}\n"
                "  state: independent\n"
                "  Ntouching: {Ntouching}\n"
                "  Nconnected: {Nconnected}\n"
                "".format(final_ncl=self.final_ncl,
                          Npx=self.Npx,
                          Tmax=self.dpeak,
                          Ntouching=len(touching),
                          Nconnected=len(connected))]
        # NOTE: PyFITS reverses the order of coordinates, therefore we
        # output 2-1-0.
        # NOTE: Coordinates in FITS start from 1 and because arrays clmask
        # and idata in this code add an one-pixel border around the original
        # data, these two shifts cancel each other out and we can
        # output ijk directly.
        str_.append("  pixels:\n")
        str_.extend(["    {ijk[2]:>3d} {ijk[1]:>3d} {ijk[0]:>3d} {dval:.12g}\n"
                     "".format(ijk=px[0], dval=px[1])
                     for px in self.pixels])
        str_.append("  touching:\n")
        str_.extend(["    {final_ncl:>3d} {dval:.12g}\n"
                     "".format(final_ncl=t[0].final_ncl, dval=t[1])
                     for t in touching])
        str_.append("  connected:\n")
        str_.extend(["    {final_ncl:>3d} {dval:.12g}\n"
                     "".format(final_ncl=c[0].final_ncl, dval=c[1])
                     for c in connected])
        str_ = "".join(str_)

        return str_


class Pixel(object):

    """Pixel within the data cube."""

    # Relative map for neighbouring pixels.
    # (Sorry for breaking PEP 8 here, it just looks better this way.) :)
    neigh_map = np.array([( 0,  0, +1),
                          ( 0,  0, -1),
                          ( 0, +1,  0),
                          ( 0, -1,  0),
                          (+1,  0,  0),
                          (-1,  0,  0)], dtype=int)

    def __init__(self, ijk, idata, clmask, clumps):
        # (i,j,k) coordinates
        self.ijk = np.array(ijk, dtype=int)
        # (x,y,z) coordinates
        self.xyz = np.array(ijk, dtype=float)
        # Input data cube
        self.idata = idata
        # Pixel mask (to which clump the pixels belong)
        self.clmask = clmask
        # List of clumps
        self.clumps = clumps
        # Data value
        self.dval = idata[tuple(self.ijk)]

    def get_neighbours(self):
        """Find neighbours touching at this pixel.

        The list of neighbours is sorted: clumps with lower ncl first.
        """
        neighbours = []
        for shift in self.neigh_map:
            ncl = self.clmask[tuple(self.ijk + shift)]
            if ncl > -1:
                neighbour = self.clumps[ncl].get_merger()
                if neighbour not in neighbours:
                    neighbours.append(neighbour)
        neighbours.sort(key=lambda clump: clump.ncl)
        return neighbours

    def get_grandparents(self, neighbours=None):
        """Find grandparents (parent's parent's... parent) among neighbours.

        The grandparents are searched for in the list of neighbours.  If the
        list is not provided, it will be taken from self.get_neighbours().

        The list of grandparents is sorted: clumps with lower ncl first.
        """
        if neighbours is None:
            neighbours = self.get_neighbours()
        grandparents = []
        for neighbour in neighbours:
            grandparent = neighbour.get_grandparent()
            if grandparent not in grandparents:
                grandparents.append(grandparent)
        grandparents.sort(key=lambda clump: clump.ncl)
        return grandparents

    def dist2(self, other):
        """Return square of the distance to other pixel or clump."""
        return ((self.xyz-other.xyz)**2).sum()

    def addto(self, clump):
        """Add pixel to clump."""
        clump.Npx += 1
        clump.pixels.append([self.ijk, self.dval])
        clump.xyz = ((self.dval*self.xyz + clump.wxyz*clump.xyz)/
                     (self.dval + clump.wxyz))
        clump.wxyz += self.dval
        clump.sumd += self.dval


if __name__ == "__main__":
    sys.exit(main())
