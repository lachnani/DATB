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

    RectRicErr_Plot(log, path, tag)
    DoeErr_Plot(log, path, tag)
    DeeErr_Plot(log, path, tag)
    RectClroeErr_Plot(log, path, tag)
    CurvClroeErr_Plot(log, path, tag)
    FilterStatus_Plot(log, path, tag)
    MeasResidual_Plot(log, path, tag)
    ChiefCov_Plot(log, path, tag)
    DeputyCov_Plot(log, path, tag)
    RelCov_Plot(log, path, tag)
    
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