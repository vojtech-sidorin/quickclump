#!/usr/bin/env python
# Unit tests for Dendroclump
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

import dclump


class TestMain(unittest.TestCase):
    """Test function main, i.e. the main Dendroclump's functionality"""

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
        """Run Dendroclump on sample files and check results."""
        for f in self.SAMPLE_FILES_PREFIXES:
            # Run main() on the sample file.
            ifits = os.path.join(self.FIXTURES_DIR, f + ".fits")
            test_ofits = os.path.join(self.TMP_DIR, f + ".clumps.fits")
            test_otext = os.path.join(self.TMP_DIR, f + ".clumps.txt")
            dclump.main("--ofits {ofits} --otext {otext} {ifits} --silent"
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


class TestParseArgs(unittest.TestCase):
    """Test function parse_args."""

    # Map command line argument strings to dicts with expected parsed
    # options.  For the test to pass, all the expected options with
    # the expected values must be found in the namespace returned by
    # function parse_args.
    ARGS_MAP = [
            ["my_fits.fits",
                {
                "ifits": "my_fits.fits",
                "dTleaf": None,
                "Tcutoff": None,
                "Npxmin": dclump.DEFAULT_NPXMIN,
                "ofits": None,
                "otext": None,
                "verbose": dclump.DEFAULT_VERBOSE
                }
            ],
            ["my_fits.fits --dTleaf 1.234 --Tcutoff 2.345 --Npxmin 3 "
                "--ofits out.fits --otext clumps.txt -vvv",
                {
                "dTleaf": 1.234,
                "Tcutoff": 2.345,
                "Npxmin": 3,
                "ofits": "out.fits",
                "otext": "clumps.txt",
                "verbose": 3
                }
            ],
            ["my_fits.fits --silent",
                {
                "ifits": "my_fits.fits",
                "verbose": dclump.SILENT_VERBOSE
                }
            ]
            ]

    def test_correct_args(self):
        for args, expected in self.ARGS_MAP:
            parsed = vars(dclump.parse_args(args.split()))
            for parameter, value in expected.items():
                self.assertIn(parameter, parsed)
                self.assertEqual(value, parsed[parameter])


class TestLoadIdata(unittest.TestCase):
    """Test function load_idata."""

    NON_3D_FITS = ["./fixtures/1d.fits",
                   "./fixtures/2d.fits",
                   "./fixtures/4d.fits"]

    THREE_DIM_FITS = ["./fixtures/3d.fits",
                      "./fixtures/rand_normal.fits",
                      "./fixtures/rand_uniform.fits"]

    def test_non_3d_fits(self):
        for filename in self.NON_3D_FITS:
            self.assertRaises(dclump.InputDataError,
                              dclump.load_idata, filename)

    def test_3d_fits(self):
        for filename in self.THREE_DIM_FITS:
            idata = dclump.load_idata(filename)
            self.assertEqual(idata.ndim, 3)

    def test_if_border_minus_inf(self):
        """Test if idata are surrounded with -inf border."""
        for filename in self.THREE_DIM_FITS:
            idata = dclump.load_idata(filename)
            self.assertTrue(np.all(np.isneginf(idata[0,:,:])))
            self.assertTrue(np.all(np.isneginf(idata[-1,:,:])))
            self.assertTrue(np.all(np.isneginf(idata[:,0,:])))
            self.assertTrue(np.all(np.isneginf(idata[:,-1,:])))
            self.assertTrue(np.all(np.isneginf(idata[:,:,0])))
            self.assertTrue(np.all(np.isneginf(idata[:,:,-1])))


class TestSetDefaults(unittest.TestCase):
    """Test function set_defaults."""

    class OptionsContainer(object):
        def __init__(self):
            self.ifits = "my_fits.fits"
            self.ofits = "output.fits"
            self.otext = "clumps.txt"
            self.dTleaf = 1.234
            self.Tcutoff = 2.345

    class InputDataContainer(object):
        def __init__(self):
            self.ndim = 3

    def test_dont_return_same_object(self):
        options = self.OptionsContainer()
        idata = self.InputDataContainer()
        self.assertIsNot(dclump.set_defaults(options, idata), options)

    def test_all_options_set(self):
        options = self.OptionsContainer()
        idata = self.InputDataContainer()
        expected = self.OptionsContainer()
        self.assertEqual(vars(dclump.set_defaults(options, idata)),
                         vars(expected))

    def test_ofits_not_set(self):
        options = self.OptionsContainer()
        options.ifits = "my_file.fits"
        options.ofits = None
        expected = self.OptionsContainer()
        expected.ifits = "my_file.fits"
        expected.ofits = "my_file.clumps.fits"
        idata = self.InputDataContainer()
        self.assertEqual(vars(dclump.set_defaults(options, idata)),
                         vars(expected))

    def test_otext_not_set(self):
        options = self.OptionsContainer()
        options.ifits = "my_file.fits"
        options.otext = None
        expected = self.OptionsContainer()
        expected.ifits = "my_file.fits"
        expected.otext = "my_file.clumps.txt"
        idata = self.InputDataContainer()
        self.assertEqual(vars(dclump.set_defaults(options, idata)),
                         vars(expected))


class TestCheckOptions(unittest.TestCase):
    """Test function check_options."""

    def setUp(self):
        class Empty(object): pass
        self.options = Empty()

    def test_missing_options(self):
        # Passing no options.
        self.assertRaises(AssertionError, dclump.check_options, self.options)
        # Missing any required option.
        required_options = ("dTleaf", "Tcutoff")
        for option in required_options:
            setattr(self.options, option, 1.)
            self.assertRaises(AssertionError, dclump.check_options,
                              self.options)
            delattr(self.options, option)

    def test_correct_values(self):
        self.options.Tcutoff = 1.
        self.options.dTleaf = 1.
        self.assertIsNone(dclump.check_options(self.options))

    def test_incorrect_values(self):
        # negative Tcutoff
        self.options.Tcutoff = -1.
        self.options.dTleaf = 1.
        self.assertRaises(dclump.OutOfBoundsError, dclump.check_options,
                          self.options)
        # negative dTleaf
        self.options.Tcutoff = 1.
        self.options.dTleaf = -1.
        self.assertRaises(dclump.OutOfBoundsError, dclump.check_options,
                          self.options)
        # negative Tcutoff and dTleaf
        self.options.Tcutoff = -1.
        self.options.dTleaf = -1.
        self.assertRaises(dclump.OutOfBoundsError, dclump.check_options,
                          self.options)
        # nan Tcutoff
        self.options.Tcutoff = float("nan")
        self.options.dTleaf = 1.
        self.assertRaises(dclump.OutOfBoundsError, dclump.check_options,
                          self.options)
        # -inf Tcutoff
        self.options.Tcutoff = float("-inf")
        self.options.dTleaf = 1.
        self.assertRaises(dclump.OutOfBoundsError, dclump.check_options,
                          self.options)


if __name__ == "__main__":
    unittest.main(verbosity=2)
