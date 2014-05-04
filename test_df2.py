#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Tests for df2
# Author: Vojtech Sidorin

import unittest

import df2

class TestCheckOptions(unittest.TestCase):
    """Test function check_options."""

    def setUp(self):
        class Empty(object): pass
        self.options = Empty()

    def test_correct_values(self):
        self.options.Tcutoff = 1.
        self.options.dTleaf = 1.
        self.assertIsNone(df2.check_options(self.options))

    def test_incorrect_values(self):
        # negative Tcutoff
        self.options.Tcutoff = -1.
        self.options.dTleaf = 1.
        self.assertRaises(df2.OutOfBoundsError, df2.check_options, self.options)
        # negative dTleaf
        self.options.Tcutoff = 1.
        self.options.dTleaf = -1.
        self.assertRaises(df2.OutOfBoundsError, df2.check_options, self.options)
        # negative Tcutoff and dTleaf
        self.options.Tcutoff = -1.
        self.options.dTleaf = -1.
        self.assertRaises(df2.OutOfBoundsError, df2.check_options, self.options)
        # nan Tcutoff
        self.options.Tcutoff = float("nan")
        self.options.dTleaf = 1.
        self.assertRaises(df2.OutOfBoundsError, df2.check_options, self.options)
        # -inf Tcutoff
        self.options.Tcutoff = float("-inf")
        self.options.dTleaf = 1.
        self.assertRaises(df2.OutOfBoundsError, df2.check_options, self.options)


if __name__ == "__main__":
    unittest.main(verbosity=2)

