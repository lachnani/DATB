# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:08:59 2025

@author: Hakim Lachnani
"""

import sys
sys.path.insert(0, '..')

import unittest
import orbit as orb
import formation as frm
import numpy as np
from kinematics import kinematicsUtils as uKin


class TestOrbit(unittest.TestCase):
    """Tests for the Orbit Class"""
    
    def setUp(self):
        """
        Create orbits for test in all modules. 
        """
        self.epoch = 0;
        self.state = np.array([6738,0.00059,np.deg2rad(51.6384),
                          np.deg2rad(144.1824),np.deg2rad(311.2672),
                          np.deg2rad(187.8960)]) # ISS on 24 February 2025 16:43:22
        self.orb1 = orb.Orbit(self.epoch,self.state,stateType = 'STATE_KEPEL')
        self.doe = np.array([5,0.0001,np.deg2rad(1),
                          np.deg2rad(0),np.deg2rad(0),
                          np.deg2rad(0.05)])
        self.frm1 = frm.Formation(self.orb1, None, self.doe, frmType='FORMATION_CHIEF_ANCHOR',
                                  relStateType='RELSTATE_DOE')
        
    def test_init(self):
        """
        Test the conversions between states which occurs at initialization

        """
        oe2 = self.frm1.deputy.oe
        relPosRectRic = np.zeros((3,))
        relVelRectRic = np.zeros((3,))
        uKin.rv2ric(self.frm1.chief.r, self.frm1.chief.v, \
            self.frm1.deputy.r, self.frm1.deputy.v, relPosRectRic, relVelRectRic)
        
        self.assertAlmostEqual(np.all(oe2-self.state), np.all(self.doe), 4)
        self.assertAlmostEqual(np.all(relPosRectRic), np.all(self.frm1.relPosRectRic), 4)
        self.assertAlmostEqual(np.all(relVelRectRic), np.all(self.frm1.relVelRectRic), 4)

        
    def test_prop(self):
        """
        Test orbit propagation

        """
        dt = 1 # 1 Hz propagation
        t = 0
        for tStep in range(1,24*3600):
            # 1 day of propagation
            t += dt
            self.frm1.propagate(dt)
            
        # Since perturbations are off, all elements should match the initial
        # set except mean anomaly. Furthermore, the relative doe state should not
        # have changed except delta mean anomaly (since da =/= 0)
        self.assertEqual(self.frm1.t, t)
        for stateIdx in range(0,4):
            self.assertAlmostEqual(self.frm1.chief.oe[stateIdx], self.state[stateIdx],4)
            self.assertAlmostEqual(self.frm1.doe[stateIdx], self.doe[stateIdx], 4) 
            
        # Chief and deputy ephemerides should match exactly
        self.assertEqual(np.all(self.frm1.chief.moon.r), np.all(self.frm1.deputy.moon.r))
        self.assertEqual(np.all(self.frm1.chief.sun.r), np.all(self.frm1.deputy.sun.r))
            
    def test_rect2curv(self):
        """
        Test rectilinear to curvilinear transformations

        """
        r = np.array([1000,0,0])
        v = np.array([0,1,0])
        relPosCurvRic0 = np.array([0,np.linalg.norm(r)*np.pi/2,0])
        relVelCurvRic0 = np.array([0,-np.linalg.norm(v)*np.pi/2,0])
        relPosRectRic = np.zeros((3,))
        relVelRectRic = np.zeros((3,))
        uKin.curvRic2rectRic(r, v, relPosCurvRic0, relVelCurvRic0, relPosRectRic, relVelRectRic)
        relPosCurvRic = np.zeros((3,))
        relVelCurvRic = np.zeros((3,))
        uKin.rectRic2curvRic(r, v, relPosRectRic, relVelRectRic, relPosCurvRic, relVelCurvRic)
        self.assertAlmostEqual(np.all(relPosCurvRic), np.all(relPosCurvRic0))
        self.assertAlmostEqual(np.all(relVelCurvRic), np.all(relVelCurvRic0))
        
        relPosCurvRic0 = np.array([0,0,np.linalg.norm(r)*np.pi/2])
        relVelCurvRic0 = np.array([0,0,0])
        uKin.curvRic2rectRic(r, v, relPosCurvRic0, relVelCurvRic0, relPosRectRic, relVelRectRic)
        uKin.rectRic2curvRic(r, v, relPosRectRic, relVelRectRic, relPosCurvRic, relVelCurvRic)
        self.assertAlmostEqual(np.all(relPosCurvRic), np.all(relPosCurvRic0))
        self.assertAlmostEqual(np.all(relVelCurvRic), np.all(relVelCurvRic0))
        
        relPosCurvRic0 = np.array([10,0,0])
        relVelCurvRic0 = np.array([0,0,0])
        uKin.curvRic2rectRic(r, v, relPosCurvRic0, relVelCurvRic0, relPosRectRic, relVelRectRic)
        uKin.rectRic2curvRic(r, v, relPosRectRic, relVelRectRic, relPosCurvRic, relVelCurvRic)
        self.assertAlmostEqual(np.all(relPosCurvRic), np.all(relPosCurvRic0))
        self.assertAlmostEqual(np.all(relVelCurvRic), np.all(relVelCurvRic0))
        
    def test_clroe(self):
        """
        Test CLROE conversions

        """
        # V-Bar Case
        relPosRectRic0 = np.array([0.0,1000.0,0.0])
        relVelRectRic0 = np.array([0.0,0.0,0.0])
        clroe0 = np.array([0.0,0.0,0.0,1000.0,0.0,0.0])
        clroe = np.zeros((6,))
        uKin.ric2clroe(relPosRectRic0, relVelRectRic0, self.frm1.chief.meanMotion, 0.0, clroe)
        self.assertAlmostEqual(np.all(clroe), np.all(clroe0))
        
        relPosRectRic = np.zeros((3,))
        relVelRectRic = np.zeros((3,))
        uKin.clroe2ric(clroe, self.frm1.chief.meanMotion, 0.0, relPosRectRic, relVelRectRic)
        self.assertAlmostEqual(np.all(relPosRectRic), np.all(relPosRectRic0))
        self.assertAlmostEqual(np.all(relVelRectRic), np.all(relVelRectRic0))
        
        # 2:1 Ellipse Case
        clroe0 = np.array([10,0.0,0.0,0.0,0.0,0.0])
        relPosRectRic = np.zeros((3,))
        relVelRectRic = np.zeros((3,))
        uKin.clroe2ric(clroe0, self.frm1.chief.meanMotion, 0.0, relPosRectRic, relVelRectRic)
        self.assertNotAlmostEqual(relPosRectRic[0], 0.0)
        self.assertAlmostEqual(relPosRectRic[1], 0.0)
        self.assertAlmostEqual(relPosRectRic[2], 0.0)
        self.assertAlmostEqual(relVelRectRic[0], 0.0)
        self.assertNotAlmostEqual(relVelRectRic[1], 0.0)
        self.assertAlmostEqual(relVelRectRic[2], 0.0)
        
        
        
if __name__ == '__main__':
    unittest.main()