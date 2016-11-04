#!/usr/bin/env python
# Tests for Quickclump.
#
# Copyright 2016 Vojtech Sidorin <vojtech.sidorin@gmail.com>
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
import shutil
import sys
import tempfile
import unittest

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

import qc


SAMPLE_FILES_DIR = "test_samples"


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
                "Npxmin": qc.DEFAULT_NPXMIN,
                "ofits": None,
                "otext": None,
                "verbose": qc.DEFAULT_VERBOSE
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
                "verbose": qc.SILENT_VERBOSE
                }
            ]
            ]

    def test_correct_args(self):
        for args, expected in self.ARGS_MAP:
            parsed = vars(qc.parse_args(args.split()))
            for parameter, value in expected.items():
                self.assertIn(parameter, parsed)
                self.assertEqual(value, parsed[parameter])


class TestLoadIdata(unittest.TestCase):
    """Test function load_idata."""

    NON_3D_FITS = [os.path.join(SAMPLE_FILES_DIR, "1d.fits"),
                   os.path.join(SAMPLE_FILES_DIR, "2d.fits"),
                   os.path.join(SAMPLE_FILES_DIR, "4d.fits")]

    THREE_DIM_FITS = [os.path.join(SAMPLE_FILES_DIR, "3d.fits"),
                      os.path.join(SAMPLE_FILES_DIR, "rand_normal.fits"),
                      os.path.join(SAMPLE_FILES_DIR, "rand_uniform.fits")]

    def test_non_3d_fits(self):
        for filename in self.NON_3D_FITS:
            self.assertRaises(qc.InputDataError, qc.load_idata, filename)

    def test_3d_fits(self):
        for filename in self.THREE_DIM_FITS:
            idata = qc.load_idata(filename)
            self.assertEqual(idata.ndim, 3)

    def test_if_border_minus_inf(self):
        """Test if idata are surrounded with -inf border."""
        for filename in self.THREE_DIM_FITS:
            idata = qc.load_idata(filename)
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
        self.assertIsNot(qc.set_defaults(options, idata), options)

    def test_all_options_set(self):
        options = self.OptionsContainer()
        idata = self.InputDataContainer()
        expected = self.OptionsContainer()
        self.assertEqual(vars(qc.set_defaults(options, idata)), vars(expected))

    def test_ofits_not_set(self):
        options = self.OptionsContainer()
        options.ifits = "my_file.fits"
        options.ofits = None
        expected = self.OptionsContainer()
        expected.ifits = "my_file.fits"
        expected.ofits = "my_file.clumps.fits"
        idata = self.InputDataContainer()
        self.assertEqual(vars(qc.set_defaults(options, idata)), vars(expected))

    def test_otext_not_set(self):
        options = self.OptionsContainer()
        options.ifits = "my_file.fits"
        options.otext = None
        expected = self.OptionsContainer()
        expected.ifits = "my_file.fits"
        expected.otext = "my_file.clumps.txt"
        idata = self.InputDataContainer()
        self.assertEqual(vars(qc.set_defaults(options, idata)), vars(expected))


class TestCheckOptions(unittest.TestCase):
    """Test function check_options."""

    def setUp(self):
        class Empty(object): pass
        self.options = Empty()

    def test_missing_options(self):
        # Passing no options.
        self.assertRaises(AssertionError, qc.check_options, self.options)
        # Missing any required option.
        required_options = ("dTleaf", "Tcutoff")
        for option in required_options:
            setattr(self.options, option, 1.)
            self.assertRaises(AssertionError, qc.check_options, self.options)
            delattr(self.options, option)

    def test_correct_values(self):
        correct_values = [1, 1., 0.1, 257, 257., 1e2, int(1e215), float("inf")]
        for self.options.Tcutoff in correct_values:
            for self.options.dTleaf in correct_values:
                self.assertIsNone(qc.check_options(self.options))

    def test_incorrect_values(self):
        correct_value = 1.
        incorrect_values = [-1, -1., -1e2, float("-inf"), float("nan")]
        # Correct Tcutoff & incorrect dTleaf.
        self.options.Tcutoff = correct_value
        for self.options.dTleaf in incorrect_values:
            self.assertRaises(qc.OutOfBoundsError, qc.check_options,
                              self.options)
        # Correct dTleaf & incorrect Tcutoff.
        self.options.dTleaf = correct_value
        for self.options.Tcutoff in incorrect_values:
            self.assertRaises(qc.OutOfBoundsError, qc.check_options,
                              self.options)
        # Incorrect both.
        for self.options.Tcutoff in incorrect_values:
            for self.options.dTleaf in incorrect_values:
                self.assertRaises(qc.OutOfBoundsError, qc.check_options,
                                  self.options)


class TestRegression(unittest.TestCase):
    """Regression test: Test if the program produces expected output."""

    # Use these sample input FITS files. (Prefixes without the .fits extension.)
    SAMPLE_IFITS_PREFIXES = ["rand_normal",
                             "rand_uniform"]

    def setUp(self):
        # Test outputs will be stored in this tmp dir.
        self.tmpd = tempfile.mkdtemp(prefix="qc-")

    def tearDown(self):
        shutil.rmtree(self.tmpd)

    def test_on_sample_input_files(self):
        """Run main() on sample files and check results."""
        for f in self.SAMPLE_IFITS_PREFIXES:

            # Run main() on the sample file.
            ifits = os.path.join(SAMPLE_FILES_DIR, f + ".fits")
            result_ofits = os.path.join(self.tmpd, f + ".clumps.fits")
            result_otext = os.path.join(self.tmpd, f + ".clumps.txt")
            qc.main("--ofits {ofits} --otext {otext} {ifits} --silent"
                    .format(ofits=result_ofits, otext=result_otext, ifits=ifits)
                    .split())

            # Compare FITS results.
            expected_ofits = os.path.join(SAMPLE_FILES_DIR, f + ".clumps.fits")
            with fits.open(expected_ofits) as g:
                expected_odata = g[0].data
                expected_header = g[0].header
            with fits.open(result_ofits) as h:
                result_odata = h[0].data
                result_header = h[0].header
            # ... data.
            self.assertEqual(list(expected_odata.flatten()),
                             list(result_odata.flatten()),
                             msg="Data in FITS files '{0}' and '{1}' differ."
                             .format(expected_ofits, result_ofits))
            # ... header.
            # Before comparison, remove the card DATE from the header.  This
            # card will contain a timestamp of the file creation and should be
            # excluded from the comparison.
            del expected_header["DATE"]
            del result_header["DATE"]
            self.assertSequenceEqual(expected_header, result_header)

            # Compare TXT results.
            expected_otext = os.path.join(SAMPLE_FILES_DIR, f + ".clumps.txt")
            with open(expected_otext, "rb") as g:
                expected_otext_contents = g.read()
            with open(result_otext, "rb") as h:
                result_otext_contents = h.read()
            self.assertEqual(expected_otext_contents, result_otext_contents,
                             msg="Files '{0}' and {1} differ."
                             .format(expected_otext, result_otext))


if __name__ == "__main__":
    unittest.main(verbosity=2)
