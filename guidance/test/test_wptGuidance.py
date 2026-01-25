# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:08:59 2025

@author: Hakim Lachnani
"""

import sys
sys.path.insert(0, '..')

import unittest
from planning import wptTbl
import linearizedModels as lm
import wptGuidance as wg
import numpy as np


class TestWptGuidance(unittest.TestCase):
    """Tests for Waypoint Guidance Functions"""
    
    def setUp(self):
        self.stm = lm.hcw(1/(24*3600)) # Intialize in GEO
        self.r0 = np.array((10,15,-8))
        self.v0 = np.zeros((3,))
        self.tf = 1000
    
    def testTwoImpulse(self):
        """
        Tests two-impulse burn calculation
        """
        # Null test. If free-flight trajectory, the burns should be null
        rf, vf = self.stm.prop(0, self.r0, self.v0, self.tf) # Natural free-flight trajectory
        dv0, dvf = wg.twoImpulseBurn(self.stm.Phi(0,self.tf), self.r0, self.v0, rf, vf)
        self.assertEqual(np.all(dv0),  np.all(np.zeros((3,))))
        self.assertEqual(np.all(dvf),  np.all(np.zeros((3,))))
        
        # Active test. Set the final position and velocity to 0 (docking)
        dv0, dvf = wg.twoImpulseBurn(self.stm.Phi(0,self.tf), self.r0, self.v0, np.zeros((3,)), np.zeros((3,)))
        v0 = self.v0 + dv0
        rf, vf = self.stm.prop(0, self.r0, v0, self.tf) # Post-burn free-flight trajectory
        vf = vf + dvf
        self.assertEqual(np.all(rf),  np.all(np.zeros((3,))))
        self.assertEqual(np.all(vf),  np.all(np.zeros((3,))))
        
    def testWptBurn(self):
        """ 
        Tests waypoint burn calculation
        """
        # Null test. If free-flight trajectory, burn should be null
        rf, vf = self.stm.prop(0, self.r0, self.v0, self.tf) # Natural free-flight trajectory
        dv0 = wg.wptBurn(self.stm.Phi(0,self.tf), self.r0, self.v0, rf)
        self.assertEqual(np.all(dv0),  np.all(np.zeros((3,))))
        
        # Active test. Set the final position to 0 (docking)
        dv0 = wg.wptBurn(self.stm.Phi(0,self.tf), self.r0, self.v0, np.zeros((3,)))
        v0 = self.v0 + dv0
        rf, vf = self.stm.prop(0, self.r0, v0, self.tf) # Post-burn free-flight trajectory
        self.assertEqual(np.all(rf),  np.all(np.zeros((3,))))
        
    def testWptBurnTable(self):
        """ 
        Tests waypoint burn table
        """
        # Initialize waypoint table
        tbl = wptTbl.waypointTable()
        tbl.genClroeWpts(np.array([4.0, 0.0, 0.0, 10.0, 2.0, 0.0]), 1/(24*3600), self.tf, 24*3600*8*np.pi, 16)
        
        # Compute burn table. First two burns should force onto natural trajectory.
        # Subsequent burns should be near zero
        burnT, burnTbl = wg.wptBurnTable(self.stm, tbl, 0, self.r0, self.v0)
        self.assertEqual(burnTbl.shape[1], tbl.numWpts + 1)
        self.assertGreater(np.linalg.norm(burnTbl[:,0]), 1e-3)
        self.assertGreater(np.linalg.norm(burnTbl[:,1]), 1e-3)
        for burn in range(2,17):
            self.assertLess(np.linalg.norm(burnTbl[:,burn]), 1e-16)
        
        
if __name__ == '__main__':
    unittest.main()