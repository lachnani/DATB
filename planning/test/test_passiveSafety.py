# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:08:59 2025

@author: Hakim Lachnani
"""

import sys
sys.path.insert(0, '..')

import unittest
import passiveSafety as ps
import numpy as np


class TestKorim(unittest.TestCase):
    """Tests for KORIM"""
    
    def test(self):
        """
        Create CLROE sets 
        A0 = L[0]
        alpha = L[1]
        xOff = L[2]
        yOff = L[3]
        B0 = L[4]
        beta = L[5]
        """
        # Case 1: Offset Radial
        self.assertEqual(ps.korim(np.array([3,0,10,0,3,0]), 2,  np.array([1,7,1])), True) # True
        self.assertEqual(ps.korim(np.array([3,0,0,0,3,0]), 2,  np.array([1,7,1])), False) # False
        
        # Case 2: Offset radial and cross-track phasing
        self.assertEqual(ps.korim(np.array([3,0,0,0,3,np.pi/2]), 2,  np.array([1,7,1])), True) # True
        self.assertEqual(ps.korim(np.array([3,0,0,0,3,0]), 2,  np.array([1,7,1])), False) # False
        
        # Case 3: Drifting Away
        self.assertEqual(ps.korim(np.array([3,0,-1,15,3,0]), 2,  np.array([1,7,1])), True) # True
        self.assertEqual(ps.korim(np.array([3,0,1,15,3,0]), 2,  np.array([1,7,1])), False) # False
        
        # Case 4: Monotonic In-Track
        self.assertEqual(ps.korim(np.array([0,0,1,5,3,0]), 2,  np.array([2,1,1])), True) # True
        self.assertEqual(ps.korim(np.array([0,0,1,7,3,0]), 2,  np.array([2,1,1])), False) # False

        # Case 5: I non-monotonic; extrema contained
        self.assertEqual(ps.korim(np.array([3,0,0,0,3,0]), 2,  np.array([1,1,1])), True) # True
        self.assertEqual(ps.korim(np.array([3,0,0,0,3,0]), 2,  np.array([1,7,1])), False) # False
        
        # Case 6: I non-monotonic; extrema not contained
        self.assertEqual(ps.korim(np.array([1,0,0,0,3,0]), 2,  np.array([1,1,1])), True) # True
        self.assertEqual(ps.korim(np.array([1,0,0,0,3,np.pi/2]), 2,  np.array([1,1,1])), False) # False
        
if __name__ == '__main__':
    unittest.main()