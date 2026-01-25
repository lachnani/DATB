# -*- coding: utf-8 -*-
"""
Created on Fri Jan 16 20:13:49 2026

@author: Hakim Lachnani
"""

import numpy as np
import scipy

class hcw():
    """ 
    Hill-Clohessy-Wiltshire Class.
    
    State vector is given as RIC position and velocities
    
    """
    
    def __init__(self, meanMotion):
        """
        Initialize Class

        """
        self.n = meanMotion
        
    def A(self):
        """
        State Matrix A from the dynamics equation xDot = A*x + B*u

        Returns
        -------
        A : 6x6 state matrix.

        """
        A = np.array((
            (0, 0, 0, 1, 0, 0),
            (0, 0, 0, 0, 1, 0),
            (0, 0, 0, 0, 0, 1),
            (3*self.n**2, 0, 0, 0, 2*self.n, 0),
            (0, 0, 0, -2*self.n, 0, 0),
            (0, 0, -self.n**2, 0, 0, 0)
            ))
        
        return A
    
    def B(self):
        """
        Input Matrix B from the dynamics equation xDot = A*x + B*u

        Returns
        -------
        B : 6x3 input matrix.

        """
        B = np.array((
            (0, 0, 0),
            (0, 0, 0),
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1)
            ))
        
        return B
        
        
    def Phi(self, t0, tf):
        """
        State Transition Matrix

        Parameters
        ----------
        t0 : initial time.
        tf : final time.

        Returns
        -------
        Phi : 6x6 STM.

        """
        dt = tf-t0
        nt = self.n*dt
        s = np.sin(nt)
        c = np.cos(nt)
        nInv = 1/self.n
        
        Phi = np.array((
            (4-3*c, 0, 0, s*nInv, 2*(1-c)*nInv, 0),
            (6*(s-nt), 1, 0, 2*(c-1)*nInv, (4*s-3*nt)*nInv, 0),
            (0, 0, c, 0, 0, s*nInv),
            (3*self.n*s, 0, 0, c, 2*s, 0),
            (6*self.n*(c-1), 0, 0, -2*s, 4*c-3, 0),
            (0, 0, -self.n*s, 0, 0, c)
            ))
        
        return Phi
    
    def prop(self, t0, r0, v0, tf):
        """
        Propagates state

        Parameters
        ----------
        t0 : initial time
        r0 : initial relative position in RIC.
        v0 : initial relative velocity in RIC.
        tf : final time

        Returns
        -------
        rf : final relative position in RIC.
        vf : final relative velocity in RIC.

        """
        
        Phi = self.Phi(t0, tf)
        rf, vf = multPhi(Phi, r0, v0)
        
        return rf, vf
    
class clroe():
    """ 
    Circular Linerized Relative Orbit Element Class.
    
    State vector is given as CLROEs state L:
        A0 = L[0]
        alpha = L[1]
        xOff = L[2]
        yOff = L[3]
        B0 = L[4]
        beta = L[5]
    
    """
    
    def __init__(self, meanMotion):
        """
        Initialize Class

        """
        self.n = meanMotion
        
    def A(self):
        """
        State Matrix A from the dynamics equation xDot = A*x + B*u

        Returns
        -------
        A : 6x6 state matrix.

        """
        return np.zeros((6,6))
    
    def B(self, L, t):
        """
        Input Matrix B from the dynamics equation xDot = A*x + B*u
        Ref: Bennett and Schaub, "Continuous-Time Modeling and Control using 
        Linearized Relative Orbit Elements," AAS/AIAA Astrodynamics Specialists
        Conference, 2015

        Returns
        -------
        B : 6x3 input matrix.

        """
        nt = self.meanMotion*t
        B = np.zeros((6,6))
        
        if np.abs(L[0]) > np.finfo(np.float64).eps:
            s = np.sin(L[1] + nt)
            c = np.cos(L[1] + nt)
            B[0,0] = -s
            B[0,1] = -2*c
            B[1,0] = -c/L[0]
            B[1,1] = 2*s/L[0]
        else:
            # Assume alpha = 0, 
            B[0,0] = -np.sin(nt)
            B[0,1] = -2*np.cos(nt)
            
        if np.abs(L[4]) > np.finfo(np.float64).eps:
            s = np.sin(L[5] + nt)
            c = np.cos(L[5] + nt)
            B[2,2] = -s
            B[3,2] = -c/L[4]
        else:
            # Assume beta = 0
            B[2,2] = -np.sin(nt)
            
        B[4,1] = 2
        B[5,0] = -2
        B[5,1] = 3*nt
        
        B = B / self.meanMotion
        
        return B
        
        
    def Phi(self):
        """
        State Transition Matrix

        Returns
        -------
        Phi : 6x6 STM.

        """        
        return np.eye(6)
    
    def prop(self, L0):
        """
        Propagates state

        Parameters
        ----------
        L0 : initial state.

        Returns
        -------
        Lf : final state.

        """        
        return L0
    
def multPhi(Phi, r0, v0):
    rf = np.matmul(Phi[0:3,0:3],r0) + np.matmul(Phi[0:3,3:6],v0)
    vf = np.matmul(Phi[3:6,0:3],r0) + + np.matmul(Phi[3:6,3:6],v0)
    
    return rf, vf

def xDot(A, x, B = np.zeros((6,3)), u = np.zeros((3,))):
    return np.matmul(A, x) + np.matmul(B, u)

def integrate(t0, tf, A, x, B = np.zeros((6,3)), u = np.zeros((3,))):
    def f(t, x, A, B, u): return xDot(A, x, B, u)
    sol = scipy.solve_ivp(f, [t0, tf], x, args = (A, B, u))
    return sol.y[:,-1]
    