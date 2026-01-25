# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:08:59 2025

@author: Hakim Lachnani
"""

import sys
sys.path.insert(0, '..')

import unittest
import wptTbl
import numpy as np


class TestWptTbl(unittest.TestCase):
    """Tests for Waypoint Table Generators"""
    
    def testInit(self):
        """
        Tests basic waypoint functionality
        """
        # Initialize Table
        tbl = wptTbl.waypointTable()
        
        # Add waypoints
        t = np.linspace(1, 3, 3)
        pos = np.array([[1,2,3],[4,5,6],[7,8,9]])
        tbl.appendWpts(t, pos)
        self.assertEqual(tbl.numWpts, 3)
        self.assertEqual(tbl.t[2], 3) 
        self.assertEqual(np.all(tbl.relPosRic[:,2]),  np.all(np.array([7,8,9]))) 
        
        # Increment
        tbl.increment()
        self.assertEqual(tbl.numWpts, 2)
        self.assertEqual(tbl.t[1], 3) 
        self.assertEqual(np.all(tbl.relPosRic[:,1]),  np.all(np.array([7,8,9]))) 
        
    def testClroeWpts(self):
        """
        Tests CLROE Waypoint Generation
        """
        # Initialize Table
        tbl = wptTbl.waypointTable()
        
        # Add waypoints
        tbl.genClroeWpts(np.array([4.0, 0.0, 0.0, 10.0, 2.0, 0.0]), 1, 0, 8*np.pi, 16)
        self.assertEqual(tbl.numWpts, 16)
        self.assertEqual(tbl.t[-1], 8*np.pi) 
        self.assertEqual(np.all(tbl.relPosRic[:,-1]),  np.all(np.array([4.0,10.0,2.0]))) 
        
    def testRptWpts(self):
        """
        Tests Rpt Waypoint Generation
        """
        # Initialize Table
        tbl = wptTbl.waypointTable()
        
        # Add waypoints
        tbl.genRptWpts(np.array([1.0,2.0,3.0]), 0, 10, 3)
        self.assertEqual(tbl.numWpts, 3)
        self.assertEqual(tbl.t[-1], 20) 
        self.assertEqual(np.all(tbl.relPosRic[:,-1]),  np.all(np.array([1.0,2.0,3.0]))) 
        
    def testFmcWpts(self):
        """
        Tests FMC Waypoint Generation
        """
        L = np.array([4.0, 0.0, 1.0, 10.0, 2.0, 0.0])
        
        # Single Rev; variable range
        tbl = wptTbl.waypointTable()
        tbl.genFmcWpts(L, 1, 0, 1, 8, numRevs=1, constantRange=False)
        self.assertEqual(tbl.numWpts, 8)
        self.assertEqual(tbl.t[-1], 2*np.pi) 
        self.assertEqual(np.all(tbl.relPosRic[:,-1]),  np.all(np.array([4.0,10.0,2.0]))) 
        
        # Multiple Revs; variable range
        tbl = wptTbl.waypointTable()
        tbl.genFmcWpts(L, 1, 0, 1, 8, numRevs=3, constantRange=False)
        self.assertEqual(tbl.numWpts, 24)
        self.assertEqual(tbl.t[-1], 6*np.pi) 
        self.assertEqual(np.all(tbl.relPosRic[:,-1]),  np.all(np.array([4.0,10.0,2.0]))) 
        
        # Multiple Revs; constant range
        tbl = wptTbl.waypointTable()
        tbl.genFmcWpts(L, 1, 0, 1, 8, numRevs=3, constantRange=True)
        self.assertEqual(tbl.numWpts, 24)
        self.assertEqual(tbl.t[-1], 6*np.pi) 
        self.assertEqual(np.all(tbl.relPosRic[:,-1]),  np.all(np.array([4.0,10.0,2.0])/np.linalg.norm(np.array([4.0,10.0,2.0])))) 
        
if __name__ == '__main__':
    unittest.main()