#!/usr/bin/env python
# Quickclump - identify clumps within a 3D FITS datacube.
#
# Copyright 2019 Vojtech Sidorin <vojtech.sidorin@gmail.com>
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

    $ python qc.py my_datacube.fits

To show usage help, type

    $ python qc.py -h

To run in the interactive mode, type

     $ python
     >>> import qc
     >>> qc.main(["my_datacube.fits"])

For more information, see README.md.
"""


# Semantic versioning; see <http://semver.org/>.
__version__ = "1.4.1"

import argparse
from collections import deque
import copy
import datetime
import logging
import os
import sys

# Import FITS IO.
# NOTE: PyFITS was merged into Astropy.
try:
    from astropy.io import fits
except ImportError:
    try:
        import pyfits as fits
    except ImportError:
        sys.exit("Error: Cannot find any supported FITS IO package.  "
                 "Have you installed 'astropy' or 'pyfits'?")
import numpy as np

# Global settings.
DEFAULT_NPXMIN = 5
DEFAULT_LOGLEVEL = "INFO"
VERBOSE_LOGLEVEL = "DEBUG"
SILENT_LOGLEVEL = "CRITICAL"
# Relative map of a pixel's neighbourhood.
# pragma pylint: disable=bad-whitespace
PIXEL_NEIGHBOURHOOD = (( 0,  0, +1),
                       ( 0,  0, -1),
                       ( 0, +1,  0),
                       ( 0, -1,  0),
                       (+1,  0,  0),
                       (-1,  0,  0))
# pragma pylint: enable=bad-whitespace

# Set a global logger.
LOGGER = logging.getLogger(__name__)
FORMATTER = logging.Formatter("%(levelname)s - %(message)s")
STDOUT_HANDLER = logging.StreamHandler(sys.stdout)
STDOUT_HANDLER.setFormatter(FORMATTER)
STDOUT_FILTER = logging.Filter()
STDOUT_FILTER.filter = lambda rec: rec.levelno < logging.ERROR
STDOUT_HANDLER.addFilter(STDOUT_FILTER)
STDERR_HANDLER = logging.StreamHandler(sys.stderr)
STDERR_HANDLER.setFormatter(FORMATTER)
STDERR_HANDLER.setLevel(logging.ERROR)
LOGGER.addHandler(STDOUT_HANDLER)
LOGGER.addHandler(STDERR_HANDLER)


def main(argv=None):
    """The main entry point."""
    try:
        _main(argv=argv)
    except (IOError, InputDataError, OutOfBoundsError) as exc:
        LOGGER.error(exc)
        return 1
    else:
        return 0


def _main(argv=None):
    # NOTE: If argv is None, the arguments from sys.argv will be used instead.
    options = parse_args(argv)

    # Load the input data (a FITS datacube).
    idata = load_idata(options.ifits)

    options = set_defaults(options, idata)
    check_options(options)

    # Initialise the clumps mask: pixels labeled with the number of the
    # clump to which they belong.
    clmask = np.empty(idata.shape, dtype="int32")
    clmask[:] = -1
    # NOTE: dtype will be reviewed later, before saving the output into a
    # FITS file, and changed to a smaller int sufficient for holding the
    # number of the found clumps.
    # NOTE: Initially, the clumps are numbered from 0 on, -1 meaning no
    # clump owning the pixel.  The final numbering, stored in clumps'
    # atribute final_ncl, will start from 1 with 0 meaning no clumps owning
    # the pixel.

    # The discovered clumps will be collected in this list.
    clumps = []

    LOGGER.debug("Finding clumps.")
    find_all_clumps(idata, clmask, clumps, options)

    LOGGER.debug("Merging small clumps.")
    merge_small_clumps(clumps, options.Npxmin)

    LOGGER.debug("Renumbering clumps.")
    final_clumps_count = renumber_clumps(clumps, options.Npxmin)
    renumber_clmask(clmask, clumps)

    # NOTE: The clumps have now set their final labels/numbers, stored in the
    # attribute final_ncl.
    # NOTE: The final_ncl of clumps that are too smal, with Npx < Npxmin, has
    # been set to 0.

    if options.ofits.strip().upper() != "NONE":
        LOGGER.debug("Writing output FITS.")
        write_ofits(options.ofits, clmask, final_clumps_count, options)

    if options.otext.strip().upper() != "NONE":
        LOGGER.debug("Writing output text file.")
        write_otext(options.otext, clumps, options)

    LOGGER.info("%i clumps found.", final_clumps_count)


def parse_args(argv=None):
    """Parse arguments."""
    parser = argparse.ArgumentParser(
        description="Identifies clumps within a 3D FITS datacube.")
    parser.add_argument("ifits", help="FITS file where to search for clumps.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--dTleaf",
        type=float,
        help="Minimal depth of a valley separating adjacent clumps.  Clumps "
             "separated by a valley that is shallower will be merged "
             "together.  Must be > 0.  (default: 3*sig_noise)"
        )
    parser.add_argument(
        "--Tcutoff",
        type=float,
        help="Minimal data value to consider.  Pixels with lower values won't "
             "be processed.  Must be > 0.  (default: 3*sig_noise)"
        )
    parser.add_argument(
        "--Npxmin",
        type=int,
        default=DEFAULT_NPXMIN,
        help="Minimal size of a clump in pixels.  Smaller clumps will be "
             "either merged to an adjacent clumps or deleted.  "
             "(default: %(default)s)"
        )
    parser.add_argument(
        "--ofits",
        help="FITS file where the found clumps will be saved.  If OFITS "
             "exists, it will be overwritten.  If set to 'None' (case doesn't "
             "matter), OFITS file won't be written.  "
             "(default: ifits with modified extension '.clumps.fits')"
        )
    parser.add_argument(
        "--otext",
        help="Text file where the found clumps will be saved in a "
             "human-readable form.  If OTEXT exists, it will be overwritten.  "
             "If set to 'None' (case insensitive), OTEXT file won't be "
             "written.  This will speed up the program's execution.  On the "
             "other hand, the OTEXT file is needed for the construction of a "
             "dendrogram.  (default: ifits with modified extension "
             "'.clumps.txt')")
    parser.add_argument(
        "--loglevel",
        dest="loglevel",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Threshold determining which types of messages will be shown. "
             "Only messages with severity greater than or equal to a given "
             "level will be shown. (default: %(default)s)",
        default=DEFAULT_LOGLEVEL
        )
    parser.add_argument(
        "--verbose", "-v",
        action="store_const",
        const=VERBOSE_LOGLEVEL,
        dest="loglevel",
        help="Increase verbosity. (Set loglevel to {0}.)"
             "".format(VERBOSE_LOGLEVEL)
        )
    parser.add_argument(
        "--silent",
        action="store_const",
        const=SILENT_LOGLEVEL,
        dest="loglevel",
        help="Decrease verbosity. (Set loglevel to {0}.)"
             "".format(SILENT_LOGLEVEL)
        )
    args = parser.parse_args(argv)
    LOGGER.setLevel(args.loglevel)
    return args


def load_idata(ifits):

    """Load and preprocess input FITS data."""

    # Load the first HDU from the FITS (HDU = header data unit).
    with fits.open(ifits) as f:  # pylint: disable=invalid-name
        # f[0] == the first HDU in the file.
        idata = f[0].data  # pylint: disable=no-member

    # If idata is of integral type, convert to float.  This is to allow
    # the addition of a border around the idata with values of -inf.
    if np.issubdtype(idata.dtype, np.integer):
        idata = idata.astype("f8")

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

    # Create a shallow copy of options.
    new_options = copy.copy(options)

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
        LOGGER.debug("Options dTleaf and/or Tcutoff was not set. "
                     "Estimating from the input data.")

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

        # Check if estimation of std_noise succeeded.
        if (not np.isfinite(std_noise)) or (std_noise <= 0.):
            raise OutOfBoundsError(
                "Estimation of std_noise from input data failed.  "
                "Got value '{0}'.  "
                "Is the input FITS data valid/reasonable?".format(std_noise))

        # Set dTleaf.
        if new_options.dTleaf is None:
            LOGGER.debug("Setting dTleaf to %f (= 3*std_noise = 3*%f)",
                         3.*std_noise, std_noise)
            new_options.dTleaf = 3.*std_noise

        # Set Tcutoff.
        if new_options.Tcutoff is None:
            LOGGER.debug("Setting Tcutoff to %f (= 3*std_noise = 3*%f)",
                         3.*std_noise, std_noise)
            new_options.Tcutoff = 3.*std_noise

    return new_options


def check_options(options):
    """Check values of dTleaf and Tcutoff."""
    assert hasattr(options, "dTleaf")
    assert hasattr(options, "Tcutoff")
    if not options.dTleaf > 0.:
        raise OutOfBoundsError("'dTleaf' must be > 0. It is {0}."
                               .format(options.dTleaf))
    if not options.Tcutoff > 0.:
        raise OutOfBoundsError("'Tcutoff' must be > 0. It is {0}."
                               .format(options.Tcutoff))


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

    # Sort flattened keys of idata array.
    sorted_px_fkeys = idata.argsort(axis=None)

    # Find clumps -- loop over sorted keys of pixels starting at maximum.
    ncl = -1  # Initialise clump index.
    assert options.Tcutoff > idata[0, 0, 0]
    for px_fkey in sorted_px_fkeys[::-1]:

        # Derive n-dim px key
        px_key = np.unravel_index(px_fkey, idata.shape)

        # Get data value
        dval = idata[px_key]

        # Skip NANs
        if dval is np.nan:
            continue

        # Terminate if dval < Tcutoff.  The keys are sorted so we don't need
        # to process the remaining pixels.
        if dval < options.Tcutoff:
            break

        # Initialize pixel
        px = Pixel(px_key, dval)

        # Find neighbours (clumps touching at this pixel)
        neighbours = sorted(px.get_neighbours(clmask, clumps),
                            key=lambda clump: clump.dpeak, reverse=True)

        if not neighbours:
            # No neighbour --> Make a new clump
            ncl += 1
            clumps.append(Clump(ncl, px))
            clmask[px_key] = ncl
        elif len(neighbours) == 1:
            # One neighbour --> Add pixel to it
            clmask[px_key] = neighbours[0].ncl
            neighbours[0].add_px(px)
        else:
            # More neighbours --> Merge/connect them
            # NOTE: There are two things to do now:
            #
            #  (1) Add the pixel to a clump.
            #  (2) Update the properties of the neighbouring clumps:
            #       (a) Merge too short clumps.
            #       (b) Update touching lists.
            #       (c) Connect root parents.

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
            clmask[px_key] = merger.ncl
            merger.add_px(px)

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

            # (2c) Connect root parents.
            # This is done so that when one side of a U-shaped clump group
            # should be rejected, e.g. because of having too little pixels,
            # that side should rather be merged to the other side than be
            # deleted.
            root_parents = sorted(px.get_root_parents(neighbours),
                                  key=lambda clump: clump.dpeak, reverse=True)
            for i in range(1, len(root_parents)):
                root_parent = root_parents[i]
                root_parent.parent = root_parents[0]
                dist2_min = root_parent.dist2(root_parent.parent)
                for j in range(1, i):
                    dist2 = root_parent.dist2(root_parents[j])
                    if dist2 < dist2_min:
                        root_parent.parent = root_parents[j]
                        dist2_min = dist2


def merge_small_clumps(clumps, Npxmin):
    """Merge clumps with too little pixels.

    Clumps with too little pixels (< Npxmin) will be merged to their
    parents.  Clumps will also be merged, if their parents have too little
    pixels.  Merging starts from the "bottom", i.e. clumps with the lowest
    dpeak values will be processed first.
    """
    for clump in reversed(clumps):
        if clump.merged:
            # Already merged --> skip
            continue
        elif clump.parent is clump:
            # Solitary/orphan clump --> skip
            continue
        elif clump.npx < Npxmin:
            # Too small clump --> merge to its parent
            clump.merge_to_parent()
        elif clump.parent.merger.npx < Npxmin:
            # Too small parent --> merge clump to it
            clump.merge_to_parent()


def renumber_clumps(clumps, Npxmin):
    """Renumber clumps taking into account mergers and Npxmin limit.

    Set clumps' final_ncl so that:
     - The numbering starts from 1.
     - Clumps that are merged or have too little pixels (<Npxmin) are
       excluded.

    Return the final count of clumps: last new_ncl.
    """
    new_ncl = 0
    for clump in clumps:
        if clump.merged or clump.npx < Npxmin:
            continue
        else:
            new_ncl += 1
            clump.final_ncl = new_ncl
    return new_ncl


def renumber_clmask(clmask, clumps):
    """Renumber clmask according to clumps' final_ncl."""
    clmask[:] = 0
    for clump in clumps:
        if clump.merged or clump.final_ncl is None:
            continue
        else:
            for px in clump.pixels:
                clmask[tuple(px.ijk)] = clump.final_ncl


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
    ohdu.header["COMMENT"] = ("File created by Quickclump {v}."
                              .format(v=__version__))
    ohdu.header["COMMENT"] = ("Original data file: '{ifits}'"
                              .format(ifits=os.path.basename(options.ifits)))
    ohdu.header["COMMENT"] = ("Quickclump was run with following parameters:")
    ohdu.header["COMMENT"] = ("  dTleaf={dTleaf:.12g}"
                              .format(dTleaf=options.dTleaf))
    ohdu.header["COMMENT"] = ("  Tcutoff={Tcutoff:.12g}"
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

    with open(otext, "w") as f:  # pylint: disable=invalid-name
        # Output the header line.
        # NOTE: The Nlevels value is set to 1000 only for the output to be
        # compatible with DENDROFIND.  It has no meaning for Quickclump.
        f.write("# Nlevels = 1000 Tcutoff = {options.Tcutoff:.12g} dTleaf = "
                "{options.dTleaf:.12g} Npxmin = {options.Npxmin}\n"
                .format(options=options))
        # Output the clumps.
        for clump in clumps:
            # Output only clumps that were not rejected (final_ncl is not None)
            # or merged.
            if (not clump.merged) and (clump.final_ncl is not None):
                f.write(str(clump))
                f.write("\n")


class InputDataError(Exception):
    """Error in the input FITS data."""


class OutOfBoundsError(Exception):
    """Error signaling that a variable is out of expected bounds."""


class PixelLike(object):
    """Base class for pixel-like objects."""
    xyz = None
    def dist2(self, other):
        """Return square of the distance to other object."""
        return ((self.xyz-other.xyz)**2).sum()


class Clump(PixelLike):

    """Clump found within the data cube."""

    def __init__(self, ncl, px):
        """Build a new clump with first pixel px.

        Positional arguments:

         ncl -- number/label for this clump
         px  -- the first pixel of the clump (peak)
        """
        super(Clump, self).__init__()
        # Clump label
        self.ncl = ncl
        # Final clump label to be set during last renumbering.
        self.final_ncl = None
        # List of pixels belonging to the clump.
        self._pixels = [px]
        # Parent to which the clump is connected (the nearest clump with a
        # higher dpeak touching this clump).
        self._parent = self
        # Whether the clump was merged to its parent.
        self._merged = False
        # Other clumps which touch this one.
        # {clump_reference: dval_at_which_they_touch}
        self.touching = {}
        # Other clumps connected with this one.
        # {clump: dval at which connect}
        self.connected = {}
        # (x,y,z) coordinates of the clump:  The weighted average with the
        # weight equal to the clump's pixels data values.  Note the xyz changes
        # as new pixels are being added to the clump.
        self.xyz = px.xyz
        # The weight of the clump's xyz coordinates.
        self.wxyz = px.dval
        # The peak data value.
        self.dpeak = px.dval

    @property
    def pixels(self):
        """Return the pixels assigned to the clump."""
        return self._pixels

    @property
    def npx(self):
        """Return the number of pixels assigned to the clump."""
        return len(self._pixels)

    @property
    def parent(self):
        """Return the clump's parent if it exists or the clump itself."""
        if self._parent is self:
            return self
        return self._parent.merger

    @parent.setter
    def parent(self, new_parent):
        self._parent = new_parent.merger

    @parent.deleter
    def parent(self):
        self._parent = self

    @property
    def root_parent(self):
        """Return self or parent's parent's... parent."""
        if self._parent is self:
            return self
        return self._parent.root_parent

    @property
    def merger(self):
        """Return self or clump to which this clump merges."""
        if self._merged:
            return self.parent
        return self

    @property
    def merged(self):
        """Return true if this clump was merged, false otherwise."""
        return self._merged

    def merge_to_parent(self):
        """Merge clump to its parent."""
        parent = self.parent
        assert self is not parent, "Attempt to merge clump to itself."
        assert not self.merged, "Attempt to merge already merged clump."
        assert not parent.merged, "Attempt to merge to merged clump."
        parent.pixels.extend(self.pixels)
        parent.xyz = ((parent.wxyz*parent.xyz + self.wxyz*self.xyz)/
                      (parent.wxyz + self.wxyz))
        parent.wxyz += self.wxyz
        for clump, touching_at_dval in self.touching.items():
            parent.update_touching(clump, touching_at_dval)
        self._merged = True

    def update_touching(self, other, dval):
        """Add clump other to the dict of touching clumps.

        Positional arguments:

         other -- The clump which touches this one.  Note that the other clump
                  is expanded to the merger.
         dval  -- The data value of the pixel which connects the two clumps.
        """
        exp_other = other.merger
        if exp_other is not self:
            if ((exp_other not in self.touching) or
                    (dval > self.touching[exp_other])):
                self.touching.update({exp_other: dval})

    def compact_touching(self):
        """Compact the touching dict.

        (1) Remove references to rejected clumps (with final_ncl == None).
        (2) Remove references to the clump itself. (That may result from
            merging a touched clump to this clump.)
        (3) Ensure the touching dict is unique (solving mergers).  If more
            references to the same clump are found, sets the touching_at_dval
            to the highest value.
        """
        new_touching = {}
        for clump, touching_at_dval in self.touching.items():
            exp_clump = clump.merger
            if exp_clump.final_ncl is None:
                continue
            elif exp_clump is self:
                continue
            elif ((exp_clump not in new_touching) or
                  (touching_at_dval > new_touching[exp_clump])):
                new_touching.update({exp_clump: touching_at_dval})
        self.touching = new_touching

    def get_connected(self):

        """Return clumps connected to this clump.

        Connected clumps are all clumps that either touch this clump directly
        or indirectly through other connected clumps.  This structure is used
        for building a dendrogram.  In other words, connected clumps make up a
        graph data structure.  We now want to find those clumps (nodes) -- i.e.
        we want to discover the whole graph.

        This method saves the last returned value in self.connected to improve
        the performance of successive graph traversals.

        Returns a dict with clumps connected to this clump in the form of
        {clump: connects_at_dval, ...}, where connects_at_dval is the data
        value at which a given clump connects.  The clump itself is not
        contained in the dict.
        """

        # Init the queue of clumps to explore; start with self.merger.
        # Queue format: [[clump, dval], ...], where dval is the data value at
        # which a given clump connects.
        # NOTE: Only mergers are expected in the queue.
        queue = deque()
        queue.append((self.merger, self.merger.dpeak))

        # Dict with connected clumps that will be returned.
        # NOTE: Only expanded clumps are expected in the dict.
        connected = dict(queue)

        # Find all connected clumps (discover the whole graph).
        while queue:
            next_clump = queue.popleft()  # Breadth-first traversal.
            focused_clump = next_clump[0]
            focused_valley = next_clump[1]
            assert not focused_clump.merged, \
                "Only expanded clumps are expected in the queue."
            if focused_clump.connected:
                # Reuse what has been already discovered.
                for child, child_valley in focused_clump.connected.items():
                    min_valley = min(focused_valley, child_valley)
                    if child in connected:
                        min_valley = max(min_valley, connected[child])
                    connected.update({child: min_valley})
            else:
                for child, child_valley in focused_clump.touching.items():
                    # Get the merger.
                    child_merger = child.merger
                    # Get the minimal data value along the path.
                    min_valley = min(focused_valley, child_valley)
                    if child_merger in connected:
                        # Rediscovered clump; update if found a better/"higher"
                        # path (with greater minimal dval along it).
                        if min_valley > connected[child_merger]:
                            queue.append((child_merger, min_valley))
                            connected[child_merger] = min_valley
                    else:
                        # Newly discovered clump
                        queue.append((child_merger, min_valley))
                        connected.update({child_merger: min_valley})

        # Exclude the clump itself from connected.
        del connected[self.merger]
        self.connected = connected
        return connected

    def add_px(self, px):
        """Add pixel to the clump."""
        self.pixels.append(px)
        self.xyz = ((px.dval*px.xyz + self.wxyz*self.xyz)/
                    (px.dval + self.wxyz))
        self.wxyz += px.dval

    def __str__(self):

        """Return textual representation of the clump.

        The output text is compatible with the original dendrofind's textual
        output and should work with the script for plotting dendrograms
        (df_dendrogram.py).
        """

        self.compact_touching()

        # Sort touching clumps (order by dval_at_which_they_touch, ncl).
        touching = list(self.touching.items())
        touching.sort(key=lambda x: (-x[1], x[0].final_ncl))

        # Get connected clumps.
        # Then sort them (by dval_at_which_they_connect, ncl).
        connected = list(self.get_connected().items())
        connected.sort(key=lambda x: (-x[1], x[0].final_ncl))

        # Sort list of pixels (order by dval, k, j, i).
        self.pixels.sort(key=lambda px:
                         (-px.dval, px.ijk[2], px.ijk[1], px.ijk[0]))

        # Generate str_ to be returned.
        str_ = ["clump: {final_ncl}\n"
                "  Npx: {Npx}\n"
                "  Tmax: {Tmax:.12g}\n"
                "  state: independent\n"
                "  Ntouching: {Ntouching}\n"
                "  Nconnected: {Nconnected}\n"
                "".format(final_ncl=self.final_ncl,
                          Npx=self.npx,
                          Tmax=float(self.dpeak),
                          Ntouching=len(touching),
                          Nconnected=len(connected))]
        # NOTE: FITS IO reverses the order of coordinates, therefore we
        # output 2-1-0.
        # NOTE: Coordinates in FITS start from 1 and because arrays clmask
        # and idata in this code add an one-pixel border around the original
        # data, these two shifts cancel each other out and we can
        # output ijk directly.
        str_.append("  pixels:\n")
        str_.extend(["    {ijk[2]:>3d} {ijk[1]:>3d} {ijk[0]:>3d} {dval:.12g}\n"
                     "".format(ijk=px.ijk, dval=float(px.dval))
                     for px in self.pixels])
        str_.append("  touching:\n")
        str_.extend(["    {final_ncl:>3d} {dval:.12g}\n"
                     "".format(final_ncl=t[0].final_ncl, dval=float(t[1]))
                     for t in touching])
        str_.append("  connected:\n")
        str_.extend(["    {final_ncl:>3d} {dval:.12g}\n"
                     "".format(final_ncl=c[0].final_ncl, dval=float(c[1]))
                     for c in connected])

        return "".join(str_)


class Pixel(PixelLike):

    """Pixel within the data cube."""

    def __init__(self, ijk, dval):
        super(Pixel, self).__init__()
        # (i,j,k) coordinates
        self.ijk = np.array(ijk, dtype=int)
        # Data value
        self.dval = dval

    @property
    def xyz(self):
        """XYZ coordinates of the pixel."""
        return self.ijk.astype(float)

    def get_neighbours(self, clmask, clumps):
        """Find neighbours touching at this pixel."""
        neighbours = set()
        for shift in PIXEL_NEIGHBOURHOOD:
            ncl = clmask[tuple(self.ijk + shift)]
            if ncl > -1:
                neighbours.add(clumps[ncl].merger)
        return neighbours

    def get_root_parents(self, neighbours):
        """Find root parents among neighbours."""
        return set([neighbour.root_parent for neighbour in neighbours])


if __name__ == "__main__":
    sys.exit(main())
