# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
    
    
# Line Labels
line_vec = ['X','Y','Z']
line_ric = ['R','I','C']
line_oe = ['a', 'e', 'i', 'RAAN', 'argP', 'M']
line_ee = ['a', 'l', 'P1', 'P2', 'Q1', 'Q2']
line_q = ['qs', 'qx', 'qy', 'qz']
line_doe = ['da', 'de', 'di', 'dRAAN', 'dargP', 'dM']
line_dee = ['da', 'dl', 'dP1', 'dP2', 'dQ1', 'dQ2']
line_roe = ['da', 'dlambda', 'dex', 'dey', 'dix', 'diy']
line_clroe = ['A0', 'alpha', 'xOff', 'yOff', 'B0', 'beta']
line_measParams = ['Rng', 'Rng-Rate', 'Az', 'El']
line_covDiag = ['XX', 'YY', 'ZZ']
line_covDiagRic = ['RR','II','CC']
line_meas = ['Az', 'El', 'Rng', 'Rng-Rate']
line_fswComp = ['Truth','FSW']

# Axis Labels
ax_posVel = ['Position [km]','Velocity [km/s]']
ax_posPos = ['Position [km]','Position [km]']
ax_oe = ['a [km]','e','i [deg]',r'$\Omega$ [deg]',r'$\omega$ [deg]','M [deg]']
ax_ee = ['a [km]','l [deg]','P1','P2','Q1','Q2']
ax_doe = [r'$\delta a$ [km]',r'$\delta e$',r'$\delta i$ [deg]',r'$\delta RAAN$ [deg]',r'$\delta\omega$ [deg]',r'$\delta M$ [deg]']
ax_dee = [r'$\delta a$ [km]',r'$\delta l [deg]$',r'$\delta P1$',r'$\delta P2$',r'$\delta Q1$',r'$\delta Q2$']
ax_clroe = [r'$A$ [km]',r'$\alpha$ [deg]',r'$x_{off}$ [km]',r'$y_{off}$ [km]',r'$B$ [km]',r'$\beta$ [deg]']
ax_posVelVar = ['Position [km^2]','Velocity [(km/s)^2]']
ax_posPosVar = ['Position [km^2]','Position [km^2]']

# Titles
title_oe = ['Semimajor-Axis','Eccentricity','Inclination','RAAN','Argument of Perigee','Mean Anomaly']
title_ee = ['Semimajor-Axis','Mean Longitude','P1','P2','Q1','Q2']
title_doe = ['Delta Semimajor-Axis','Delta Eccentricity','Delta Inclination','Delta RAAN','Delta Argument of Perigee','Delta Mean Anomaly']
title_dee = ['Delta Semimajor-Axis','Delta Mean Longitude','Delta P1','Delta P2','Delta Q1','Delta Q2']
title_clroe = ['In-Plane Ellipse Size','In-Plane Ellipse Phase','Radial Offset','In-Track Offset','Cross-Track Magnitude','Cross-Track Phase']   
    
    
def plotThreeByTwo(time, data, axLabels, titles, fullFigPath, suptitle = None):
    """
    3x2 plot

    """
    
    fig = plt.figure()
    if suptitle is not None:
        plt.suptitle(suptitle)
        
    for ii in range(6):
        ax = plt.subplot(3,2,ii+1)
        ax.plot(time, data[ii,:])[0]
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(axLabels[ii])
        plt.title(titles[ii])
        plt.grid()
       
    fig.tight_layout()
    
    plt.savefig(fullFigPath)
    
def plotTwoByTwo(time, data, axLabels, titles, fullFigPath, suptitle = None):
    """
    2x2 plot

    """
    
    fig = plt.figure()
    if suptitle is not None:
        plt.suptitle(suptitle)
        
    for ii in range(4):
        ax = plt.subplot(2,2,ii+1)
        ax.plot(time, data[ii,:])[0]
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(axLabels[ii])
        plt.title(titles[ii])
        plt.grid()
       
    fig.tight_layout()
    
    plt.savefig(fullFigPath)
    
def plotThreeByTwoComp(time, data, lineLabels, axLabels, titles, fullFigPath, suptitle = None):
    """
    3x2 comparison plot

    """
    
    fig = plt.figure()
    if suptitle is not None:
        plt.suptitle(suptitle)
        
    nData = len(data)
        
    for ii in range(6):
        ax = plt.subplot(3,2,ii+1)
        for jj in range(nData):
            ax.plot(time, data[jj][ii,:], label=lineLabels[jj])[0]
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(axLabels[ii])
        ax.legend()
        plt.title(titles[ii])
        plt.grid()
       
    fig.tight_layout()
    
    plt.savefig(fullFigPath)
    
def plotStack(time, data, lineLabels, axLabels, titles, fullFigPath, suptitle = None):
    """
    Stacked plot

    """
    
    fig = plt.figure()
    if suptitle is not None:
        plt.suptitle(suptitle)
        
    nStack = len(data)
        
    for iStack in range(nStack):
        ax = plt.subplot(nStack,1,iStack+1)
        nData = len(data[iStack])
        for ii in range(nData):
            ax.plot(time, data[iStack][ii,:], label=lineLabels[iStack][ii])[0]
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(axLabels[iStack])
        ax.legend()
        plt.title(titles[iStack])
        plt.grid()
       
    fig.tight_layout()
    
    plt.savefig(fullFigPath)
    