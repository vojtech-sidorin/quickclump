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

import unittest

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


if __name__ == "__main__":
    unittest.main(verbosity=2)
