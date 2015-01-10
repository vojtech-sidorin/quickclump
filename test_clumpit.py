#!/usr/bin/env python
# Unit tests for clumpit
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

import os
import sys
import shutil
import unittest
import hashlib

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

import clumpit

class TestCheckOptions(unittest.TestCase):
    """Test function check_options."""

    def setUp(self):
        class Empty(object): pass
        self.options = Empty()

    def test_missing_options(self):
        # Passing no options.
        self.assertRaises(AssertionError, clumpit.check_options, self.options)
        # Missing any required option.
        required_options = ("dTleaf", "Tcutoff")
        for option in required_options:
            setattr(self.options, option, 1.)
            self.assertRaises(AssertionError, clumpit.check_options,
                              self.options)
            delattr(self.options, option)

    def test_correct_values(self):
        self.options.Tcutoff = 1.
        self.options.dTleaf = 1.
        self.assertIsNone(clumpit.check_options(self.options))

    def test_incorrect_values(self):
        # negative Tcutoff
        self.options.Tcutoff = -1.
        self.options.dTleaf = 1.
        self.assertRaises(clumpit.OutOfBoundsError, clumpit.check_options,
                          self.options)
        # negative dTleaf
        self.options.Tcutoff = 1.
        self.options.dTleaf = -1.
        self.assertRaises(clumpit.OutOfBoundsError, clumpit.check_options,
                          self.options)
        # negative Tcutoff and dTleaf
        self.options.Tcutoff = -1.
        self.options.dTleaf = -1.
        self.assertRaises(clumpit.OutOfBoundsError, clumpit.check_options,
                          self.options)
        # nan Tcutoff
        self.options.Tcutoff = float("nan")
        self.options.dTleaf = 1.
        self.assertRaises(clumpit.OutOfBoundsError, clumpit.check_options,
                          self.options)
        # -inf Tcutoff
        self.options.Tcutoff = float("-inf")
        self.options.dTleaf = 1.
        self.assertRaises(clumpit.OutOfBoundsError, clumpit.check_options,
                          self.options)

class TestMain(unittest.TestCase):
    """Test function main, i.e. the main clumpit's functionality"""

    TMP_DIR = "./test_tmp"
    FIXTURES_DIR = "./fixtures"
    SAMPLE_FILES_PREFIXES = ["rand_normal",
                             "rand_uniform"]

    def setUp(self):
        # Make tmp directory for test outputs.  If the directory exists mkdir
        # will fail, the test will fail with an exception, and the directory
        # won't be removed by tearDown.  I.e. it won't delete the directory if
        # it already exists.
        os.mkdir(self.TMP_DIR)

    def tearDown(self):
        shutil.rmtree("test_tmp")

    def test_on_sample_input_files(self):
        """Run clumpit on sample files and check results."""
        for f in self.SAMPLE_FILES_PREFIXES:
            # Run main() on the sample file.
            ifits = os.path.join(self.FIXTURES_DIR, f + ".fits")
            test_ofits = os.path.join(self.TMP_DIR, f + ".clumps.fits")
            test_otext = os.path.join(self.TMP_DIR, f + ".clumps.txt")
            clumpit.main("--ofits {ofits} --otext {otext} {ifits} --silent"
                         .format(ofits=test_ofits, otext=test_otext,
                                 ifits=ifits)
                         .split())
            # Compare FITS results.  Compare only data, not FITS header.
            sample_ofits = os.path.join(self.FIXTURES_DIR, f + ".clumps.fits")
            with fits.open(sample_ofits) as g:
                sample_odata = g[0].data  # First HDU.
            with fits.open(test_ofits) as h:
                test_odata = h[0].data  # First HDU.
            self.assertEqual(hashlib.sha512(sample_odata).hexdigest(),
                             hashlib.sha512(test_odata).hexdigest(),
                             msg="Data in FITS files '{0}' and '{1}' differ."
                                 .format(sample_ofits, test_ofits))
            # Compare TXT results.
            sample_otext = os.path.join(self.FIXTURES_DIR, f + ".clumps.txt")
            with open(sample_otext, "rb") as g:
                sample_otext_contents = g.read()
            with open(test_otext, "rb") as h:
                test_otext_contents = h.read()
            self.assertEqual(hashlib.sha512(sample_otext_contents).hexdigest(),
                             hashlib.sha512(test_otext_contents).hexdigest(),
                             msg="Files '{0}' and {1} differ."
                                 .format(sample_otext, test_otext))


if __name__ == "__main__":
    unittest.main(verbosity=2)
