#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script by Vojtěch Sidorin.
# This script is part of the DENDROFIND package.

"""Tests df2.py."""

import df2
import unittest

class TestCheckOptions(unittest.TestCase):
    """Tests function check_options."""
    
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
        self.assertRaises(df2.Error, df2.check_options, self.options)
        

if __name__ == "__main__":
    unittest.main()
