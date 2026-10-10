# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

def visAll(log, path, tag, settings):
    """
    Plot all available plots and animations

    """
    plotAll(log, path, tag, settings)
    
def plotAll(log, path, tag, settings):
    """
    Plot all available plots

    """
    # Absolute State Plots
    ChiefOe_Plot(log, path, tag)
    DeputyOe_Plot(log, path, tag)
    # Relative State Plots
    RectRicProj_Plot(log, path, tag)
    CurvRicProj_Plot(log, path, tag)
    # Error Plots
    RectRicErr_Plot(log, path, tag)
    DoeErr_Plot(log, path, tag)
    DeeErr_Plot(log, path, tag)
    RectClroeErr_Plot(log, path, tag)
    CurvClroeErr_Plot(log, path, tag)
    # Filter Status Plot
    FilterStatus_Plot(log, path, tag)
    # Measument Residual Plot
    MeasResidual_Plot(log, path, tag)
    # Covariance Plots
    ChiefCov_Plot(log, path, tag)
    DeputyCov_Plot(log, path, tag)
    RelCov_Plot(log, path, tag)
    SmaVariance_Plot(log, path, tag)
    
def ChiefOe_Plot(log, path, tag):
    """
    Chief Keplerian Elements

    """
    
    fig_oe_plt = plt.figure()
    plt.suptitle('Chief Orbit Elements')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.oeRso.loc['a',:], label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], log.fswNavOeRso.loc['a',:], label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    plt.title('Semimajor-Axis')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], log.oeRso.loc['e',:], label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], log.fswNavOeRso.loc['e',:], label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    plt.title('Eccentricity')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['i',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeRso.loc['i',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    plt.title('Inclination')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['RAAN',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeRso.loc['RAAN',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("RAAN [deg]")
    plt.title('RAAN')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['argP',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeRso.loc['argP',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    plt.title('Argument of Perigee')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['M',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeRso.loc['M',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    plt.title('Mean Anomaly')
    ax.legend()
    plt.grid()
       
    fig_oe_plt.tight_layout()
    
    fullFigPath = path + r"\chiefOe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyOe_Plot(log, path, tag):
    """
    Deputy Keplerian Elements

    """
    
    fig_oe_plt = plt.figure()
    plt.suptitle('Deputy Orbit Elements')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.oeVeh.loc['a',:], label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], log.fswNavOeVeh.loc['a',:], label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    plt.title('Semimajor-Axis')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], log.oeVeh.loc['e',:], label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], log.fswNavOeVeh.loc['e',:], label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    plt.title('Eccentricity')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeVeh.loc['i',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeVeh.loc['i',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    plt.title('Inclination')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeVeh.loc['RAAN',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeVeh.loc['RAAN',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("RAAN [deg]")
    plt.title('RAAN')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeVeh.loc['argP',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeVeh.loc['argP',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    plt.title('Argument of Perigee')
    ax.legend()
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeVeh.loc['M',:]), label='Truth')[0]
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavOeVeh.loc['M',:]), label='FSW')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    plt.title('Mean Anomaly')
    ax.legend()
    plt.grid()
       
    fig_oe_plt.tight_layout()
    
    fullFigPath = path + r"\chiefOe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicProj_Plot(log, path, tag):
    """
    2D Rectilinear RIC plots

    """
    
    fig_rectRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    ax1.plot(log.relPosRectRic.loc['I',:], log.relPosRectRic.loc['R',:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosRectRic.loc['I',log.i], log.relPosRectRic.loc['R',log.i], color='r')
    ax1.plot(log.relPosRectRic.loc['I',:], log.fswNavRelPosRectRic.loc['R',:], linestyle='--', color='k', label='FSW')[0]
    ax1.scatter(log.relPosRectRic.loc['I',log.i], log.fswNavRelPosRectRic.loc['R',log.i], color='k')
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    ax2.plot(log.relPosRectRic.loc['I',:], log.relPosRectRic.loc['C',:], color='r', label='Truth')[0]
    ax2.scatter(log.relPosRectRic.loc['I',log.i], log.relPosRectRic.loc['C',log.i], color='r')
    ax2.plot(log.relPosRectRic.loc['I',:], log.fswNavRelPosRectRic.loc['C',:], linestyle='--', color='k', label='FSW')[0]
    ax2.scatter(log.relPosRectRic.loc['I',log.i], log.fswNavRelPosRectRic.loc['C',log.i], color='k')
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    ax3.plot(log.relPosRectRic.loc['R',:], log.relPosRectRic.loc['C',:], color='r', label='Truth')[0]
    ax3.scatter(log.relPosRectRic.loc['R',log.i], log.relPosRectRic.loc['C',log.i], color='r')
    ax3.plot(log.relPosRectRic.loc['R',:], log.fswNavRelPosRectRic.loc['C',:], linestyle='--', color='k', label='FSW')[0]
    ax3.scatter(log.relPosRectRic.loc['R',log.i], log.fswNavRelPosRectRic.loc['C',log.i], color='k')
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Cross-Track [km]")
    ax3.set_aspect('equal', adjustable='datalim')
    ax3.legend()
    plt.grid()
    
    fig_rectRicProj_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicProj_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicProj_Plot(log, path, tag):
    """
    2D Curvilinear RIC plots

    """
    
    fig_curvRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    ax1.plot(log.relPosCurvRic.loc['I',:], log.relPosCurvRic.loc['R',:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosCurvRic.loc['I',log.i], log.relPosCurvRic.loc['R',log.i], color='r')
    ax1.plot(log.relPosCurvRic.loc['I',:], log.fswNavRelPosCurvRic.loc['R',:], linestyle='--', color='k', label='FSW')[0]
    ax1.scatter(log.relPosCurvRic.loc['I',log.i], log.fswNavRelPosCurvRic.loc['R',log.i], color='k')
    ax1.set_xlabel("Curvilinear In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    ax2.plot(log.relPosCurvRic.loc['I',:], log.relPosCurvRic.loc['C',:], color='r', label='Truth')[0]
    ax2.scatter(log.relPosCurvRic.loc['I',log.i], log.relPosCurvRic.loc['C',log.i], color='r')
    ax2.plot(log.relPosCurvRic.loc['I',:], log.fswNavRelPosCurvRic.loc['C',:], linestyle='--', color='k', label='FSW')[0]
    ax2.scatter(log.relPosCurvRic.loc['I',log.i], log.fswNavRelPosCurvRic.loc['C',log.i], color='k')
    ax2.set_xlabel("Curvilinear In-Track [km]")
    ax2.set_ylabel("Curvilinear Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    ax3.plot(log.relPosCurvRic.loc['R',:], log.relPosCurvRic.loc['C',:], color='r', label='Truth')[0]
    ax3.scatter(log.relPosCurvRic.loc['R',log.i], log.relPosCurvRic.loc['C',log.i], color='r')
    ax3.plot(log.relPosCurvRic.loc['R',:], log.fswNavRelPosCurvRic.loc['C',:], linestyle='--', color='k', label='FSW')[0]
    ax3.scatter(log.relPosCurvRic.loc['R',log.i], log.fswNavRelPosCurvRic.loc['C',log.i], color='k')
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Curvilinear Cross-Track [km]")
    ax3.set_aspect('equal', adjustable='datalim')
    ax3.legend()
    plt.grid()
    
    fig_curvRicProj_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicProj_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicErr_Plot(log, path, tag):
    """
    Rectilinear RIC position and velocity errors

    """
    
    fig_rectRic_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelPosRectRicErr.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavRelPosRectRicErr.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelPosRectRicErr.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Position Error')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelVelRectRicErr.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavRelVelRectRicErr.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelVelRectRicErr.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Velocity Error')
    plt.grid()
       
    fig_rectRic_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRic_Plot(log, path, tag):
    """
    Curvilinear RIC position and velocities

    """
    
    fig_curvRic_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelPosCurvRicErr.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavRelPosCurvRicErr.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelPosCurvRicErr.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Position Error')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelVelCurvRicErr.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavRelVelCurvRicErr.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavRelVelCurvRicErr.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Velocity Error')
    plt.grid()
       
    fig_curvRic_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    
def DoeErr_Plot(log, path, tag):
    """
    Delta Keplerian Orbit element error plot

    """
    
    fig_doe_plt = plt.figure()
    plt.suptitle('Differential Keplerian Orbit Element Errors')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.fswNavDoeErr.loc['da',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta a$ [km]")
    plt.title('Delta Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], log.fswNavDoeErr.loc['de',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta e$")
    plt.title('Delta Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavDoeErr.loc['di',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta i$ [deg]")
    plt.title('Delta Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavDoeErr.loc['dRAAN',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta RAAN$ [deg]")
    plt.title('Delta RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavDoeErr.loc['dargP',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta\omega$ [deg]")
    plt.title('Delta Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavDoeErr.loc['dM',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta M$ [deg]")
    plt.title('Delta Mean Anomaly')
    plt.grid()
       
    fig_doe_plt.tight_layout()
    
    fullFigPath = path + r"\doeErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeeErr_Plot(log, path, tag):
    """
    Differential Equinoctial Element error plot

    """
    
    fig_dee_plt = plt.figure()
    plt.suptitle('Differential Equinoctial Element Errors')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.fswNavDeeErr.loc['da',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta a [km]$")
    plt.title('Delta Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavDeeErr.loc['dl',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta l [deg]$")
    plt.title('Delta Mean Longitude')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.fswNavDeeErr.loc['dP1',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta P1$")
    plt.title('Delta P1')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.fswNavDeeErr.loc['dP2',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta P2$")
    plt.title('Delta P2')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.fswNavDeeErr.loc['dQ1',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta Q1$")
    plt.title('Delta Q1')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], log.fswNavDeeErr.loc['dQ2',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta Q2$")
    plt.title('Delta Q2')
    plt.grid()
       
    fig_dee_plt.tight_layout()
    
    fullFigPath = path + r"\deeErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectClroeErr_Plot(log, path, tag):
    """
    Rectilinear CLROE Error plot

    """
    
    fig_rectClroe_plt = plt.figure()
    plt.suptitle('Rectilinear CLROE Errors')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.fswNavRectClroeErr.loc['A0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$A_0$ [km]")
    plt.title('In-Plane Ellipse Size')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavRectClroeErr.loc['alpha',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\alpha$ [deg]")
    plt.title('In-Plane Ellipse Phase')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.fswNavRectClroeErr.loc['xOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$x_{off}$ [km]")
    plt.title('Radial Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.fswNavRectClroeErr.loc['yOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$y_{off}$ [km]")
    plt.title('In-Track Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.fswNavRectClroeErr.loc['B0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$B_0$ [km]")
    plt.title('Cross-Track Magnitude')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavRectClroeErr.loc['beta',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\beta$ [deg]")
    plt.title('Cross-Track Phase')
    plt.grid()
       
    fig_rectClroe_plt.tight_layout()
    
    fullFigPath = path + r"\rectClroeErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvClroeErr_Plot(log, path, tag):
    """
    Curvilinear CLROE error plot

    """
    
    fig_curvClroe_plt = plt.figure()
    plt.suptitle('Curvilinear CLROE Errors')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.fswNavCurvClroeErr.loc['A0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$A_0$ [km]")
    plt.title('In-Plane Ellipse Size')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavCurvClroeErr.loc['alpha',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\alpha$ [deg]")
    plt.title('In-Plane Ellipse Phase')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.fswNavCurvClroeErr.loc['xOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$x_{off}$ [km]")
    plt.title('Radial Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.fswNavCurvClroeErr.loc['yOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$y_{off}$ [km]")
    plt.title('In-Track Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.fswNavCurvClroeErr.loc['B0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$B_0$ [km]")
    plt.title('Cross-Track Magnitude')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavCurvClroeErr.loc['beta',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\beta$ [deg]")
    plt.title('Cross-Track Phase')
    plt.grid()
       
    fig_curvClroe_plt.tight_layout()
    
    fullFigPath = path + r"\curvClroeErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def FilterStatus_Plot(log, path, tag):
    """
    Filter status plots

    """
    
    fig_fltrStatus_plt = plt.figure()
    plt.suptitle('Filter Statuses')
    
    ax1 = plt.subplot(5,1,1)
    ax1.plot(log.time.loc['simTime',:], log.fswNavFltrInit.loc['Status',:])[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Initialized")
    plt.grid()
    
    ax2 = plt.subplot(5,1,2)
    ax2.plot(log.time.loc['simTime',:], log.fswNavFltrConverged.loc['Status',:])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Converged")
    plt.grid()
    
    ax3 = plt.subplot(5,1,3)
    ax3.plot(log.time.loc['simTime',:], log.fswNavFltrDiverged.loc['Status',:])[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Diverged")
    plt.grid()
    
    ax4 = plt.subplot(5,1,4)
    ax4.plot(log.time.loc['simTime',:], log.fswNavFltrCorrupted.loc['Status',:])[0]
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Corrupted")
    plt.grid()
    
    ax5 = plt.subplot(5,1,5)
    ax5.plot(log.time.loc['simTime',:], log.fswNavFltrConsistent.loc['Status',:])[0]
    ax5.set_xlabel("Time [s]")
    ax5.set_ylabel("Consisten")
    plt.grid()
       
    fig_fltrStatus_plt.tight_layout()
    
    fullFigPath = path + r"\fltrStatus_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    
def MeasResidual_Plot(log, path, tag):
    """
    Measurement Residual Plots

    """
    
    fig_measResidual_plt = plt.figure()
    
    ax1 = plt.subplot(3,1,1)
    az = ax1.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavMeasRes.loc['az',:]), label='Az')[0]
    el = ax1.plot(log.time.loc['simTime',:], np.rad2deg(log.fswNavMeasRes.loc['el',:]), label='El')[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("LOS Angles [deg]")
    ax1.legend()
    plt.title('Line of Sight Angles')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    ax2.plot(log.time.loc['simTime',:], log.fswNavMeasRes.loc['rng',:])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Range [km]")
    plt.title('Range')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    ax3.plot(log.time.loc['simTime',:], log.fswNavMeasRes.loc['rngRate',:])[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Range Rate [km/s]")
    plt.title('Range Rate')
    plt.grid()
       
    fig_measResidual_plt.tight_layout()
    
    fullFigPath = path + r"\measResidual_" + tag + r".png"
    plt.savefig(fullFigPath)
    

def ChiefCov_Plot(log, path, tag):
    """
    Chief covariance diagonal elements

    """
    
    fig_chiefCov_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRsoEciDiag.loc['X',:], label='XX')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavCovRsoEciDiag.loc['Y',:], label='YY')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRsoEciDiag.loc['Z',:], label='ZZ')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km^2]")
    ax.legend()
    plt.title('Chief Inertial Position Covariance')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRsoEciDiag.loc['Xv',:], label='XX')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavCovRsoEciDiag.loc['Yv',:], label='YY')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRsoEciDiag.loc['Zv',:], label='ZZ')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [(km/s)^2]")
    ax.legend()
    plt.title('Chief Inertial Velocity Covariance')
    plt.grid()
       
    fig_chiefCov_plt.tight_layout()
    
    fullFigPath = path + r"\chiefCovariance_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyCov_Plot(log, path, tag):
    """
    Deputy covariance diagonal elements

    """
    
    fig_deputyCov_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovVehEciDiag.loc['X',:], label='XX')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavCovVehEciDiag.loc['Y',:], label='YY')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovVehEciDiag.loc['Z',:], label='ZZ')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km^2]")
    ax.legend()
    plt.title('Deputy Inertial Position Covariance')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovVehEciDiag.loc['Xv',:], label='XX')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavCovVehEciDiag.loc['Yv',:], label='YY')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovVehEciDiag.loc['Zv',:], label='ZZ')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [(km/s)^2]")
    ax.legend()
    plt.title('Deputy Inertial Velocity Covariance')
    plt.grid()
       
    fig_deputyCov_plt.tight_layout()
    
    fullFigPath = path + r"\deputyCovariance_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RelCov_Plot(log, path, tag):
    """
    Relative covariance diagonal elements

    """
    
    fig_relCov_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRelRicDiag.loc['R',:], label='RR')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavCovRelRicDiag.loc['I',:], label='II')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRelRicDiag.loc['C',:], label='CC')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km^2]")
    ax.legend()
    plt.title('Relative RIC Position Covariance')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRelRicDiag.loc['Rv',:], label='RR')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.fswNavCovRelRicDiag.loc['Iv',:], label='II')[0]
    X = ax.plot(log.time.loc['simTime',:], log.fswNavCovRelRicDiag.loc['Cv',:], label='CC')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [(km/s)^2]")
    ax.legend()
    plt.title('Realtive RIC Velocity Covariance')
    plt.grid()
       
    fig_relCov_plt.tight_layout()
    
    fullFigPath = path + r"\relativeCovariance_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def SmaVariance_Plot(log, path, tag):
    """
    Semimajor axis variance

    """
    
    fig_smaVar_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    ax.plot(log.time.loc['simTime',:], log.fswNavRsoSmaVar.loc['var',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a Variance [km^2]")
    plt.title('Chief Semimajor Axis Variance')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(log.time.loc['simTime',:], log.fswNavVehSmaVar.loc['var',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a Variance [km^2]")
    plt.title('Deputy Semimajor Axis Variance')
    plt.grid()
       
    fig_smaVar_plt.tight_layout()
    
    fullFigPath = path + r"\smaVariance_" + tag + r".png"
    plt.savefig(fullFigPath)