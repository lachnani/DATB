# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

from plotting import plotGeneric as pg 

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
    # Filter Consistency Metrics Plots
    Consistency_Plot(log, path, tag)
    
def ChiefOe_Plot(log, path, tag):
    """
    Chief Keplerian Elements

    """
    pg.plotThreeByTwoComp(
        log.simTime, 
        [log.oeRso,log.fswNavOeRso], 
        pg.line_fswComp, 
        pg.ax_oe, 
        pg.title_oe, 
        path + r"\chiefOe_" + tag + r".png",
        suptitle='Chief Orbit Elements')
    
def DeputyOe_Plot(log, path, tag):
    """
    Deputy Keplerian Elements

    """
    pg.plotThreeByTwoComp(
        log.simTime, 
        [log.oeVeh,log.fswNavOeVeh], 
        pg.line_fswComp, 
        pg.ax_oe, 
        pg.title_oe, 
        path + r"\cdeputyOe_" + tag + r".png",
        suptitle='Deputy Orbit Elements')
    
def RectRicProj_Plot(log, path, tag):
    """
    2D Rectilinear RIC plots

    """
    
    fig_rectRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    ax1.plot(log.relPosRectRic[1,:], log.relPosRectRic[0,:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosRectRic[1,log.i], log.relPosRectRic[0,log.i], color='r')
    ax1.plot(log.relPosRectRic[1,:], log.fswNavRelPosRectRic[0,:], linestyle='--', color='k', label='FSW')[0]
    ax1.scatter(log.relPosRectRic[1,log.i], log.fswNavRelPosRectRic[0,log.i], color='k')
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    ax2.plot(log.relPosRectRic[1,:], log.relPosRectRic[2,:], color='r', label='Truth')[0]
    ax2.scatter(log.relPosRectRic[1,log.i], log.relPosRectRic[2,log.i], color='r')
    ax2.plot(log.relPosRectRic[1,:], log.fswNavRelPosRectRic[2,:], linestyle='--', color='k', label='FSW')[0]
    ax2.scatter(log.relPosRectRic[1,log.i], log.fswNavRelPosRectRic[2,log.i], color='k')
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    ax3.plot(log.relPosRectRic[0,:], log.relPosRectRic[2,:], color='r', label='Truth')[0]
    ax3.scatter(log.relPosRectRic[0,log.i], log.relPosRectRic[2,log.i], color='r')
    ax3.plot(log.relPosRectRic[0,:], log.fswNavRelPosRectRic[2,:], linestyle='--', color='k', label='FSW')[0]
    ax3.scatter(log.relPosRectRic[0,log.i], log.fswNavRelPosRectRic[2,log.i], color='k')
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
    ax1.plot(log.relPosCurvRic[1,:], log.relPosCurvRic[0,:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosCurvRic[1,log.i], log.relPosCurvRic[0,log.i], color='r')
    ax1.plot(log.relPosCurvRic[1,:], log.fswNavRelPosCurvRic[0,:], linestyle='--', color='k', label='FSW')[0]
    ax1.scatter(log.relPosCurvRic[1,log.i], log.fswNavRelPosCurvRic[0,log.i], color='k')
    ax1.set_xlabel("Curvilinear In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    ax2.plot(log.relPosCurvRic[1,:], log.relPosCurvRic[2,:], color='r', label='Truth')[0]
    ax2.scatter(log.relPosCurvRic[1,log.i], log.relPosCurvRic[2,log.i], color='r')
    ax2.plot(log.relPosCurvRic[1,:], log.fswNavRelPosCurvRic[2,:], linestyle='--', color='k', label='FSW')[0]
    ax2.scatter(log.relPosCurvRic[1,log.i], log.fswNavRelPosCurvRic[2,log.i], color='k')
    ax2.set_xlabel("Curvilinear In-Track [km]")
    ax2.set_ylabel("Curvilinear Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    ax3.plot(log.relPosCurvRic[0,:], log.relPosCurvRic[2,:], color='r', label='Truth')[0]
    ax3.scatter(log.relPosCurvRic[0,log.i], log.relPosCurvRic[2,log.i], color='r')
    ax3.plot(log.relPosCurvRic[0,:], log.fswNavRelPosCurvRic[2,:], linestyle='--', color='k', label='FSW')[0]
    ax3.scatter(log.relPosCurvRic[0,log.i], log.fswNavRelPosCurvRic[2,log.i], color='k')
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
    pg.plotStack(
        log.simTime, 
        [log.fswNavRelPosRectRicErr,log.fswNavRelVelRectRicErr], 
        [pg.line_ric,pg.line_ric], 
        pg.ax_posVel, 
        ['Relative Rectilinear RIC Position Error','Relative Rectilinear RIC Velocity Error'], 
        path + r"\rectRicErr_" + tag + r".png")
    
def CurvRic_Plot(log, path, tag):
    """
    Curvilinear RIC position and velocities

    """
    pg.plotStack(
        log.simTime, 
        [log.fswNavRelPosCurvRicErr,log.fswNavRelVelCurvRicErr], 
        [pg.line_ric,pg.line_ric], 
        pg.ax_posVel, 
        ['Relative Curvilinear RIC Position Error','Relative Curvilinear RIC Velocity Error'], 
        path + r"\curvRicErr_" + tag + r".png")    
    
def DoeErr_Plot(log, path, tag):
    """
    Delta Keplerian Orbit element error plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.fswNavDoeErr, 
        pg.ax_doe, 
        pg.title_doe, 
        path + r"\doeErr_" + tag + r".png",
        suptitle='Differential Keplerian Orbit Element Errors')
    
def DeeErr_Plot(log, path, tag):
    """
    Differential Equinoctial Element error plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.fswNavDeeErr, 
        pg.ax_dee, 
        pg.title_dee, 
        path + r"\deeErr_" + tag + r".png",
        suptitle='Differential Equinoctial Element Errors')
    
def RectClroeErr_Plot(log, path, tag):
    """
    Rectilinear CLROE Error plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.fswNavRectClroeErr, 
        pg.ax_clroe, 
        pg.title_clroe, 
        path + r"\rectClroeErr_" + tag + r".png",
        suptitle='Rectilinear CLROE Errors')
    
def CurvClroeErr_Plot(log, path, tag):
    """
    Curvilinear CLROE error plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.fswNavCurvClroeErr, 
        pg.ax_clroe, 
        pg.title_clroe, 
        path + r"\curvClroeErr_" + tag + r".png",
        suptitle='Curvilinear CLROE Errors')
    
def FilterStatus_Plot(log, path, tag):
    """
    Filter status plots

    """
    
    fig_fltrStatus_plt = plt.figure()
    plt.suptitle('Filter Statuses')
    
    ax1 = plt.subplot(3,2,1)
    ax1.plot(log.simTime, log.fswNavFltrInit)[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Initialized")
    plt.grid()
    
    ax2 = plt.subplot(3,2,2)
    ax2.plot(log.simTime, log.fswNavFltrConverged)[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Converged")
    plt.grid()
    
    ax3 = plt.subplot(3,2,3)
    ax3.plot(log.simTime, log.fswNavFltrDiverged)[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Diverged")
    plt.grid()
    
    ax4 = plt.subplot(3,2,4)
    ax4.plot(log.simTime, log.fswNavFltrCorrupted)[0]
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Corrupted")
    plt.grid()
    
    ax5 = plt.subplot(3,2,5)
    ax5.plot(log.simTime, log.fswNavFltrConsistent)[0]
    ax5.set_xlabel("Time [s]")
    ax5.set_ylabel("Consistent")
    plt.grid()
    
    ax6 = plt.subplot(3,2,6)
    ax6.plot(log.simTime, log.fswNavFltrUnderweight)[0]
    ax6.set_xlabel("Time [s]")
    ax6.set_ylabel("Underweighting")
    plt.grid()
       
    fig_fltrStatus_plt.tight_layout()
    
    fullFigPath = path + r"\fltrStatus_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    
def MeasResidual_Plot(log, path, tag):
    """
    Measurement Residual Plots

    """
    
    fig_measResidual_plt = plt.figure()
    plt.suptitle('Measurement Residuals')
    
    ax1 = plt.subplot(3,1,1)
    az = ax1.plot(log.simTime, np.rad2deg(log.fswNavMeasRes[0,:]), label='Az')[0]
    el = ax1.plot(log.simTime, np.rad2deg(log.fswNavMeasRes[1,:]), label='El')[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("LOS Angles [deg]")
    ax1.legend()
    plt.title('Line of Sight Angles')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    ax2.plot(log.simTime, log.fswNavMeasRes[2,:])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Range [km]")
    plt.title('Range')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    ax3.plot(log.simTime, log.fswNavMeasRes[3,:])[0]
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
    pg.plotStack(
        log.simTime, 
        [log.fswNavCovRsoEciDiag[0:3,:],log.fswNavCovRsoEciDiag[3:6,:]], 
        [pg.line_covDiag,pg.line_covDiag], 
        pg.ax_posVelVar, 
        ['Chief Inertial Position Covariance','Chief Inertial Velocity Covariance'], 
        path + r"\chiefCovariance_" + tag + r".png") 
    
def DeputyCov_Plot(log, path, tag):
    """
    Deputy covariance diagonal elements

    """
    pg.plotStack(
        log.simTime, 
        [log.fswNavCovVehEciDiag[0:3,:],log.fswNavCovVehEciDiag[3:6,:]], 
        [pg.line_covDiag,pg.line_covDiag], 
        pg.ax_posVelVar, 
        ['Deputy Inertial Position Covariance','Deputy Inertial Velocity Covariance'], 
        path + r"\deputyCovariance_" + tag + r".png") 
    
def RelCov_Plot(log, path, tag):
    """
    Relative covariance diagonal elements

    """
    pg.plotStack(
        log.simTime, 
        [log.fswNavCovRelRicDiag[0:3,:],log.fswNavCovRelRicDiag[3:6,:]], 
        [pg.line_covDiagRic,pg.line_covDiagRic], 
        pg.ax_posVelVar, 
        ['Relative RIC Position Covariance','Realtive RIC Velocity Covariance'], 
        path + r"\relativeCovariance_" + tag + r".png") 
    
def SmaVariance_Plot(log, path, tag):
    """
    Semimajor axis variance

    """
    
    fig_smaVar_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    ax.plot(log.simTime, log.fswNavRsoSmaVar)[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a Variance [km^2]")
    plt.title('Chief Semimajor Axis Variance')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(log.simTime, log.fswNavVehSmaVar)[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a Variance [km^2]")
    plt.title('Deputy Semimajor Axis Variance')
    plt.grid()
       
    fig_smaVar_plt.tight_layout()
    
    fullFigPath = path + r"\smaVariance_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def Consistency_Plot(log, path, tag):
    """
    Filter Consistency Plots

    """
    
    fig_consistency_plt = plt.figure()
    
    ax1 = plt.subplot(2,1,1)
    ax1.plot(log.simTime, log.fswNavMeasNrmInvnSqr)[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("NIS")
    plt.title('Normalized Innovations Squared')
    plt.grid()
    
    ax2 = plt.subplot(2,1,2)
    ax2.plot(log.simTime, log.fswNavRsoNrmEstErrSqr, label='Chief')[0]
    ax2.plot(log.simTime, log.fswNavVehNrmEstErrSqr, label='Deputy')[0]
    ax2.plot(log.simTime, log.fswNavRelNrmEstErrSqr, label='Relative')[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("NEES")
    ax2.legend()
    plt.title('Normalized Estimation Error Squared')
    plt.grid()
       
    fig_consistency_plt.tight_layout()
    
    fullFigPath = path + r"\consistency_" + tag + r".png"
    plt.savefig(fullFigPath)