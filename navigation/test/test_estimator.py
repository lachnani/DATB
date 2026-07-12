# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 14:08:59 2025

@author: Hakim Lachnani
"""

import sys
sys.path.insert(0, '..')

import unittest
import numpy as np
from numpy import random as rand
import estimator as est
from dynamics import orbit as orb
from dynamics import formation
from kinematics import kinematicsUtils as uKin
import measurements as meas

class TestEstimator(unittest.TestCase):
    """Tests for the EKF Class"""
    
    def setUp(self):
        """
        Create EKF class 
        """
        self.t = 0 
        def Q(dt):
            return dt*np.diag(np.array([0,0.2,0,0.2]))
        self.Q = Q
        self.x0 = np.array([50,1,50,0])
        self.P0 = 10 * self.Q(1)
        def F(dt,x):
            A = np.eye(4)
            A[0,1] = dt 
            A[2,3] = dt 
            return A 
        self.F = F 
        def f(dt,x,u,param):
            return np.matmul(F(dt,x),x)
        self.f = f
        def z(x):
            return np.array([np.sqrt(x[0]**2+x[2]**2),np.arctan2(x[2],x[0])])
        self.z = z 
        def H(x):
            r = np.sqrt(x[0]**2+x[2]**2)
            theta = np.arctan2(x[2],x[0])
            H = np.zeros((2,4))
            H[0,0] = np.cos(theta)
            H[0,2] = np.sin(theta)
            H[1,0] = -np.sin(theta)/r
            H[1,2] = np.cos(theta)/r
            return H 
        self.H = H 
        self.R = np.diag(np.array([0.5**2,0.005**2]))
        self.S = np.zeros((2,2))
        self.ekf = est.ExtendedKalmanFilter(
            self.t, self.x0, self.P0, self.Q, self.f, self.F, self.S)
            
        
    def test_ekf(self):
        """
        Test state propagation and measurement update
        Based on: https://www.cs.cmu.edu/~16385/s17/Slides/16.4_Extended_Kalman_Filter.pdf

        """
        ### Propagation
        dt = 1
        u = np.zeros((2,))
        self.x = self.x0 + rand.multivariate_normal(np.zeros(4,),self.Q(dt))
        for ii in range(10):
            self.t = self.t + dt
            self.x = self.f(dt, self.x, u, None) + rand.multivariate_normal(np.zeros(4,),self.Q(dt))
            self.ekf.propagate(dt, u)
            
        self.assertEqual(self.t, self.ekf.t)
        self.assertEqual(self.ekf.x[0], 60)
        self.assertEqual(self.ekf.x[1], 1)
        self.assertEqual(self.ekf.x[2], 50)
        self.assertEqual(self.ekf.x[3], 0)
        self.assertGreater(np.all(np.diag(self.ekf.P)), np.all(np.diag(self.P0)))
        
        ### Measurement update
        Pminus = self.ekf.P 
        errminus = self.ekf.x - self.x
        z = self.z(self.x) + rand.multivariate_normal(np.zeros(2,),self.R)
        zHat = self.z(self.ekf.x)
        H = self.H(self.ekf.x)
        nu = z - zHat
        self.ekf.update(nu, H, self.R)
        Pplus = self.ekf.P 
        errplus = self.ekf.x - self.x
        
        self.assertEqual(self.t, self.ekf.t)   
        self.assertLess(np.abs(errplus[0]), np.abs(errminus[0]))
        self.assertLess(np.abs(errplus[2]), np.abs(errminus[2]))
        self.assertLess(Pplus[0,0], Pminus[0,0])
        self.assertLess(Pplus[2,2], Pminus[2,2])
        
    def test_diekf(self):
        """
        Test Dual Inertial EKF
        Based on: 
        Woffinden, David Charles, "Angles-Only Navigation for Autonomous 
        Orbital Rendezvous" (2008). All Graduate Theses and Dissertations. 12.
        https://digitalcommons.usu.edu/etd/12

        """
        ### Initial conditions
        tJ2000 = 0
        rc = np.zeros((3,))
        vc = np.zeros((3,))
        oec = np.array([42000.,0.,0.002,0.,0.,0.])
        meanMotion = np.sqrt(orb.MU_EARTH/oec[0]**3)
        uKin.oe2rv(orb.MU_EARTH, oec, rc, vc)
        P0 = np.block([
            [0.030**2*np.eye(3),np.zeros((3,3))],
            [np.zeros((3,3)),3.6e-6**2*np.eye(3)]])
        clroe = np.array([0.,0.,0.01,10.,0.,0.]) # Drifting co-elliptic
        relPosRic = np.zeros((3,))
        relVelRic = np.zeros((3,))
        uKin.clroe2ric(clroe, meanMotion, 0, relPosRic, relVelRic)
        rd, vd = formation.ric2rv(rc, vc, relPosRic, relVelRic)
        procVar = 0.06e-6
        dvVar = 3e-6
        measCov = (np.array([1e-3,1e-3,1e-3,1e-2])**2)*np.eye(4)
        
        ### Create formation class
        chief = orb.Orbit(tJ2000, oec, stateType = "STATE_KEPEL", pert = None, settings = None)
        frm = formation.Formation(chief, None, clroe, frmType = "FORMATION_CHIEF_ANCHOR",
                                  relStateType = "RELSTATE_RECT_CLROE", pert = None, settings = None)
        
        ### Initialize DIEKF class
        nav = est.DualInertialEKF(
            tJ2000, rc, vc, P0, rd, vd, P0, procVar, dvVar, measCov)
        
        ### Propagate through 2 hours
        tf = 2*3600
        dt = 10
        while nav.tJ2000 < tf:
            frm.propagate(dt)
            nav.propagate(dt, np.zeros((3,)))
            nav.sync()
            
        uKin.rv2oe(orb.MU_EARTH, nav.chiefPosInr, nav.chiefVelInr, oec)
        uKin.ric2clroe(nav.relPosRectRic, nav.relVelRectRic, meanMotion, 0, clroe)
        P1 = nav.ekf.P
        
        # Time is synched
        self.assertEqual(nav.tJ2000, tf)
        # Chief orbit matches truth
        self.assertAlmostEqual(np.all(oec),np.all(frm.chief.oe))
        # Relative state has not changed appreciably
        self.assertAlmostEqual(np.all(clroe),np.all(frm.rectClroe))
        # Covariance is block diagonal
        self.assertTrue(np.all(P1[6:12,0:6] == np.zeros((6,6)))) 
        self.assertTrue(np.all(P1[0:6,6:12] == np.zeros((6,6)))) 
        # Covariance has increasedj
        self.assertTrue(np.all(P1[0:6,0:6] >= P0))
        self.assertTrue(np.all(P1[6:12,6:12] >= P0))
        
        ### Ingest a measurement
        frm.dcmInr2Los = nav.dcmInr2Los
        frm.dcmRic2Los = np.matmul(frm.dcmInr2Los,np.transpose(frm.dcmInr2Ric))
        frm.az, frm.el = meas.calcAzEl(frm.chief.r, frm.deputy.r, frm.dcmInr2Los)
        nav.update(meas.get(frm, measCov), "anglesRange")
        nav.sync()
        
        P2 = nav.ekf.P
        
        # Time is synched
        self.assertEqual(nav.tJ2000, tf)
        # Covariance is no longer block diagonal
        self.assertFalse(np.all(P2[6:12,0:6] == np.zeros((6,6)))) 
        self.assertFalse(np.all(P2[0:6,6:12] == np.zeros((6,6)))) 
        # Covariance has decreased
        self.assertTrue(np.all(np.diag(P2[0:6,0:6]) <= np.diag(P1[0:6,0:6])))
        self.assertTrue(np.all(np.diag(P2[6:12,6:12]) <= np.diag(P1[0:6,0:6])))
     
    def test_rekf(self):
        """
        Test Relative EKF

        """
        ### Initial conditions
        tJ2000 = 0
        rc = np.zeros((3,))
        vc = np.zeros((3,))
        oec = np.array([42000.,0.,0.002,0.,0.,0.])
        meanMotion = np.sqrt(orb.MU_EARTH/oec[0]**3)
        uKin.oe2rv(orb.MU_EARTH, oec, rc, vc)
        P0 = np.block([
            [0.030**2*np.eye(3),np.zeros((3,3))],
            [np.zeros((3,3)),3.6e-6**2*np.eye(3)]])
        clroe = np.array([0.,0.,0.01,10.,0.,0.]) # Drifting co-elliptic
        relPosRic = np.zeros((3,))
        relVelRic = np.zeros((3,))
        uKin.clroe2ric(clroe, meanMotion, 0, relPosRic, relVelRic)
        rd, vd = formation.ric2rv(rc, vc, relPosRic, relVelRic)
        procVar = 0.06e-6
        dvVar = 3e-6
        measCov = (np.array([1e-3,1e-3,1e-3,1e-2])**2)*np.eye(4)
        
        ### Create formation class
        pert = {
            "jnum": 6,
            "solarGrav": True,
            "lunarGrav": True,
            "SRP": False,
            "drag": False,
            "Cd": 0.0,
            "normalizedArea": 0.0
            }
        chief = orb.Orbit(tJ2000, oec, stateType = "STATE_KEPEL", pert = pert, settings = None)
        frm = formation.Formation(chief, None, clroe, frmType = "FORMATION_CHIEF_ANCHOR",
                                  relStateType = "RELSTATE_RECT_CLROE", pert = pert, settings = None)
        
        ### Initialize DIEKF class
        rcErr = np.array([0.01,0.01,0.01])
        nav = est.RelativeEKF(
            tJ2000, rc+rcErr, vc, P0, rd, vd, P0, procVar, dvVar, measCov)
        
        ### Propagate through 2 hours
        tf = 2*3600
        dt = 10
        while nav.tJ2000 < tf:
            frm.propagate(dt)
            nav.propagate(dt, np.zeros((3,)))
            nav.sync()
            
        uKin.rv2oe(orb.MU_EARTH, nav.chiefPosInr, nav.chiefVelInr, oec)
        uKin.ric2clroe(nav.relPosRectRic, nav.relVelRectRic, meanMotion, 0, clroe)
        P1 = nav.ekf.P
        errPos1 = np.linalg.norm(nav.relPosRectRic - frm.relPosRectRic)
        errVel1 = np.linalg.norm(nav.relVelRectRic - frm.relVelRectRic)
        
        # Time is synched
        self.assertEqual(nav.tJ2000, tf)
        # Chief orbit matches truth
        self.assertAlmostEqual(np.all(oec),np.all(frm.chief.oe))
        # Relative state has not changed appreciably
        self.assertAlmostEqual(np.all(clroe),np.all(frm.rectClroe))
        # Covariance is block diagonal
        self.assertTrue(np.all(P1[6:12,0:6] == np.zeros((6,6)))) 
        self.assertTrue(np.all(P1[0:6,6:12] == np.zeros((6,6)))) 
        # Covariance has increasedj
        self.assertTrue(np.all(P1[0:6,0:6] >= P0))
        self.assertTrue(np.all(P1[6:12,6:12] >= P0))
        
        
        ### Ingest measurements
        frm.dcmInr2Los = nav.dcmInr2Los
        frm.dcmRic2Los = np.matmul(frm.dcmInr2Los,np.transpose(frm.dcmInr2Ric))
        frm.az, frm.el = meas.calcAzEl(frm.chief.r, frm.deputy.r, frm.dcmInr2Los)
        nav.update(meas.get(frm, np.zeros((4,4))), "anglesRangeRR")
        nav.sync()
        
        P2 = nav.ekf.P
        errPos2 = np.linalg.norm(nav.relPosRectRic - frm.relPosRectRic)
        errVel2 = np.linalg.norm(nav.relVelRectRic - frm.relVelRectRic)
        
        # Time is synched
        self.assertEqual(nav.tJ2000, tf)
        # Covariance is still block diagonal
        self.assertTrue(np.all(P2[6:12,0:6] == np.zeros((6,6)))) 
        self.assertTrue(np.all(P2[0:6,6:12] == np.zeros((6,6)))) 
        # Only relative Covariance has decreased
        self.assertTrue(np.all(P2[0:6,0:6] == P1[0:6,0:6]))
        self.assertTrue(np.all(np.diag(P2[6:12,6:12]) <= np.diag(P1[0:6,0:6])))
        # Error has decreased
        self.assertTrue(np.all(np.abs(errPos2) < np.abs(errPos1)))
        
        
if __name__ == '__main__':
    unittest.main()