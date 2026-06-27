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
from numpy import linalg as la
import estimator as est


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
        def f(dt,x,u):
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
            self.x = self.f(dt, self.x, u) + rand.multivariate_normal(np.zeros(4,),self.Q(dt))
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
        
        
if __name__ == '__main__':
    unittest.main()