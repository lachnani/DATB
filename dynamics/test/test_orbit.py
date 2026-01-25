# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:08:59 2025

@author: Hakim Lachnani
"""

import sys
sys.path.insert(0, '..')

import unittest
import orbit as orb
import numpy as np


class TestOrbit(unittest.TestCase):
    """Tests for the Orbit Class"""
    
    def setUp(self):
        """
        Create orbits for test in all modules. 
        """
        self.epoch = 0;
        self.orbPeriod = 5504
        sma = (((self.orbPeriod/(2*np.pi))**2)*orb.MU_EARTH)**(1/3)
        self.state = np.array([sma,0.00059,np.deg2rad(51.6384),
                          np.deg2rad(144.1824),np.deg2rad(311.2672),
                          np.deg2rad(187.8960)]) # ISS on 24 February 2025 16:43:22
        self.orb1 = orb.Orbit(self.epoch,self.state,stateType = "STATE_KEPEL")
        self.r1 = self.orb1.r.copy()
        self.v1 = self.orb1.v.copy()
        
    def test_init(self):
        """
        Test the conversions between states which occurs at initialization

        """
        rv1 = np.concatenate((self.orb1.r,self.orb1.v))
        ee1 = self.orb1.ee
        
        orb2 = orb.Orbit(self.epoch,rv1,stateType = "STATE_POSVEL")
        orb3 = orb.Orbit(self.epoch,ee1,stateType = "STATE_EQEL")
        
        self.assertAlmostEqual(np.all(rv1), np.all(np.concatenate((orb3.r, orb3.v))), 4)
        self.assertAlmostEqual(np.all(orb2.oe), np.all(self.orb1.oe), 4)
        self.assertAlmostEqual(np.all(ee1), np.all(orb2.ee), 4)
        self.assertAlmostEqual(np.all(orb3.oe), np.all(self.orb1.oe), 4)
        
    def test_prop(self):
        """
        Test orbit propagation

        """
        dt = 1 # 1 Hz propagation
        t = 0
        tf = 605440 # 110 orbit periods; approx 7 days
        for tStep in range(tf):
            # 1 day of propagation
            t += dt
            self.orb1.propagate(dt)
            
        # Since perturbations are off, all elements should match the initial
        # set except mean anomaly
        self.assertEqual(self.orb1.t, t)
        for stateIdx in range(0,4):
            self.assertAlmostEqual(self.orb1.oe[stateIdx], self.state[stateIdx],4)
            
        # Convert both states to pos/vel and check
        posDiff = self.orb1.r - self.r1
        maxErr = np.max(np.abs(posDiff))
        self.assertLess(maxErr, 1e-6)
        print(f"Maximum position error after 7 days of LEO propagation: {maxErr*1000:.5E} meters")
        
    def test_jPert(self):
        """
        Test J perturbations

        """
        oe1 = self.orb1.oe.copy()
        
        self.orb1.pert["jnum"] = 6
        
        dt = 1 # 1 Hz propagation
        for tStep in range(self.orbPeriod*50):
            # 1 orbit of propagation
            self.orb1.propagate(dt)
            
        # Since J perturbations are on, dominated by J2, we should see a 
        # secular change in RAAN and argP. We expect RAAN to decrease and argP 
        # to increase
        self.assertLess(self.orb1.oe[3], oe1[3])
        self.assertGreater(self.orb1.oe[4], oe1[4])
        
        self.orb1.pert["jnum"] = 0
        
    def test_ephPert(self):
        """
        Test Sun and Moon ephemeris perturbations

        """     
        self.orb1.pert["solarGrav"] = True
        self.orb1.pert["lunarGrav"] = True
        
        dt = 1 # 1 Hz propagation
        for tStep in range(self.orbPeriod*10):
            # 1 orbit of propagation
            self.orb1.propagate(dt)
            
        # No particular test here, we just want to make sure the perturbations 
        # actually doing something
        posDiff = self.orb1.r - self.r1
        maxErr = np.max(np.abs(posDiff))
        self.assertGreater(maxErr, 1e-6)
        
        self.orb1.pert["solarGrav"] = False
        self.orb1.pert["lunarGrav"] = False
        
    def test_dragPert(self):
        """
        Test drag perturbation

        """
        oe1 = self.orb1.oe.copy()
        
        self.orb1.pert["drag"] = True
        self.orb1.Cd = 2.4
        self.orb1.normA = 1/300
        
        dt = 1 # 1 Hz propagation
        for tStep in range(self.orbPeriod*50):
            # 1 orbit of propagation
            self.orb1.propagate(dt)
            
        # With drag, we expect semimajor axis to decrease on average
        self.assertLess(self.orb1.oe[0], oe1[0])
        
        self.orb1.pert["drag"] = False
        self.orb1.Cd = 0.0
        self.orb1.normA = 0.0
            
        
        
if __name__ == '__main__':
    unittest.main()