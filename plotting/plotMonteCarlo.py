# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

def visAll(mcLog, mcStats, path, tag, settings = None):
    """
    Plot all available plots and animations

    """
    plotAll(mcLog, mcStats, path, tag, settings)
    
def plotAll(mcLog, mcStats, path, tag, settings = None):
    """
    Plot all available plots

    """
    # Navigation overlay plots
    RectRicProj_Plot(mcLog, path, tag)
    CurvRicProj_Plot(mcLog, path, tag)
    
    # Error overlay plots
    RectRicPosErr_Plot(mcLog, path, tag)
    RectRicVelErr_Plot(mcLog, path, tag)
    CurvRicPosErr_Plot(mcLog, path, tag)
    CurvRicVelErr_Plot(mcLog, path, tag)
    
    # Error statistics plots
    RectRicStats_Plot(mcStats, path, tag)
    CurvRicStats_Plot(mcStats, path, tag)
    DoeStats_Plot(mcStats, path, tag)
    ChiefInrStats_Plot(mcStats, path, tag)
    ChiefOeStats_Plot(mcStats, path, tag)
    DeputyInrStats_Plot(mcStats, path, tag)
    DeputyOeStats_Plot(mcStats, path, tag)
    ConsistencyStats_Plot(mcStats, path, tag)
    
def plotComp(mcStats, mcStr, path, tag):
    """
    Plots comparison plots
    
    """
    
    RectRicStatsComp_Plot(mcStats, mcStr, path, tag)
    CurvRicStatsComp_Plot(mcStats, mcStr, path, tag)
    DoeStatsComp_Plot(mcStats, mcStr, path, tag)
    ChiefInrStatsComp_Plot(mcStats, mcStr, path, tag)
    ChiefOeStatsComp_Plot(mcStats, mcStr, path, tag)
    DeputyInrStatsComp_Plot(mcStats, mcStr, path, tag)
    DeputyOeStatsComp_Plot(mcStats, mcStr, path, tag)
    ConsistencyStatsComp_Plot(mcStats, mcStr, path, tag)
    
def RectRicProj_Plot(mcLog, path, tag):
    """
    2D Rectilinear RIC plots

    """
    
    log = mcLog[0]
    
    fig_rectRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    ax1.plot(log.relPosRectRic.loc['I',:], log.relPosRectRic.loc['R',:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosRectRic.loc['I',log.i], log.relPosRectRic.loc['R',log.i], color='r')
    for i in range(len(mcLog)):
        ax1.plot(mcLog[i].relPosRectRic.loc['I',:], mcLog[i].fswNavRelPosRectRic.loc['R',:], lw=0.5, color='k')[0]
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
    for i in range(len(mcLog)):
        ax2.plot(mcLog[i].relPosRectRic.loc['I',:], mcLog[i].fswNavRelPosRectRic.loc['C',:], lw=0.5, color='k')[0]
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
    for i in range(len(mcLog)):
        ax3.plot(mcLog[i].relPosRectRic.loc['R',:], mcLog[i].fswNavRelPosRectRic.loc['C',:], lw=0.5, color='k')[0]
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Cross-Track [km]")
    ax3.set_aspect('equal', adjustable='datalim')
    ax3.legend()
    plt.grid()
    
    fig_rectRicProj_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicProj_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicProj_Plot(mcLog, path, tag):
    """
    2D Curvilinear RIC plots

    """
    
    log = mcLog[0]
    
    fig_curvRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    ax1.plot(log.relPosCurvRic.loc['I',:], log.relPosCurvRic.loc['R',:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosCurvRic.loc['I',log.i], log.relPosCurvRic.loc['R',log.i], color='r')
    for i in range(len(mcLog)):
        ax1.plot(mcLog[i].relPosCurvRic.loc['I',:], mcLog[i].fswNavRelPosCurvRic.loc['R',:], lw=0.5, color='k')[0]
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    ax2.plot(log.relPosCurvRic.loc['I',:], log.relPosCurvRic.loc['C',:], color='r', label='Truth')[0]
    ax2.scatter(log.relPosCurvRic.loc['I',log.i], log.relPosCurvRic.loc['C',log.i], color='r')
    for i in range(len(mcLog)):
        ax2.plot(mcLog[i].relPosCurvRic.loc['I',:], mcLog[i].fswNavRelPosCurvRic.loc['C',:], lw=0.5, color='k')[0]
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    ax3.plot(log.relPosCurvRic.loc['R',:], log.relPosCurvRic.loc['C',:], color='r', label='Truth')[0]
    ax3.scatter(log.relPosCurvRic.loc['R',log.i], log.relPosCurvRic.loc['C',log.i], color='r')
    for i in range(len(mcLog)):
        ax3.plot(mcLog[i].relPosCurvRic.loc['R',:], mcLog[i].fswNavRelPosCurvRic.loc['C',:], lw=0.5, color='k')[0]
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Cross-Track [km]")
    ax3.set_aspect('equal', adjustable='datalim')
    ax3.legend()
    plt.grid()
    
    fig_curvRicProj_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicProj_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicPosErr_Plot(mcLog, path, tag):
    """
    Rectilinear RIC position errors
    
    """
    
    fig_rectRicPos_plt = plt.figure()
    plt.suptitle('Rectilinear RIC Position Errors')
    
    ax = plt.subplot(3,1,1)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelPosRectRicErr.loc['R',:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelPosRectRicErr.loc['I',:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelPosRectRicErr.loc['C',:], lw=0.5, color='C2')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Cross-Track [km]")
    plt.grid()
       
    fig_rectRicPos_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicPosErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicVelErr_Plot(mcLog, path, tag):
    """
    Rectilinear RIC velocity errors
    
    """
    
    fig_rectRicVel_plt = plt.figure()
    plt.suptitle('Rectilinear RIC Velcoity Errors')
    
    ax = plt.subplot(3,1,1)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelVelRectRicErr.loc['R',:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelVelRectRicErr.loc['I',:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelVelRectRicErr.loc['C',:], lw=0.5, color='C2')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Cross-Track [km/s]")
    plt.grid()
       
    fig_rectRicVel_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicVelErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicPosErr_Plot(mcLog, path, tag):
    """
    Curvilinear RIC position errors
    
    """
    
    fig_curvRicPos_plt = plt.figure()
    plt.suptitle('Curvilinear RIC Position Errors')
    
    ax = plt.subplot(3,1,1)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelPosCurvRicErr.loc['R',:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelPosCurvRicErr.loc['I',:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelPosCurvRicErr.loc['C',:], lw=0.5, color='C2')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Cross-Track [km]")
    plt.grid()
       
    fig_curvRicPos_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicPosErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicVelErr_Plot(mcLog, path, tag):
    """
    Curvilinear RIC velocity errors
    
    """
    
    fig_curvRicVel_plt = plt.figure()
    plt.suptitle('Curvilinear RIC Velcoity Errors')
    
    ax = plt.subplot(3,1,1)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelVelCurvRicErr.loc['R',:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelVelCurvRicErr.loc['I',:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].time.loc['simTime',:], mcLog[i].fswNavRelVelCurvRicErr.loc['C',:], lw=0.5, color='C2')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Cross-Track [km/s]")
    plt.grid()
       
    fig_curvRicVel_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicVelErr_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicStats_Plot(mcStats, path, tag):
    """
    Rectilinear RIC position and velocity statistics

    """
    
    fig_rectRic_plt = plt.figure()
    plt.suptitle(r"Rectilinear RIC 3$\sigma$ Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], 3*mcStats['relPosRectRicErrStd'][0,:], label='R')[0]
    ax.plot(mcStats['time'], 3*mcStats['relPosRectRicErrStd'][1,:], label='I')[0]
    ax.plot(mcStats['time'], 3*mcStats['relPosRectRicErrStd'][2,:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Position Error')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], 3*mcStats['relVelRectRicErrStd'][0,:], label='R')[0]
    ax.plot(mcStats['time'], 3*mcStats['relVelRectRicErrStd'][1,:], label='I')[0]
    ax.plot(mcStats['time'], 3*mcStats['relVelRectRicErrStd'][2,:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Velocity Error')
    plt.grid()
       
    fig_rectRic_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    fig_rectRic_plt = plt.figure()
    plt.suptitle(r"Rectilinear RIC Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], mcStats['relPosRectRicRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    plt.title('Relative Rectilinear RIC Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], mcStats['relVelRectRicRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    plt.title('Relative Rectilinear RIC Velocity RMSE')
    plt.grid()
       
    fig_rectRic_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicRMSE_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicStats_Plot(mcStats, path, tag):
    """
    Curvilinear RIC position and velocity statistics

    """
    
    fig_curvRic_plt = plt.figure()
    plt.suptitle(r"Curvilinear RIC 3$\sigma$ Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], 3*mcStats['relPosCurvRicErrStd'][0,:], label='R')[0]
    ax.plot(mcStats['time'], 3*mcStats['relPosCurvRicErrStd'][1,:], label='I')[0]
    ax.plot(mcStats['time'], 3*mcStats['relPosCurvRicErrStd'][2,:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Position Error')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], 3*mcStats['relVelCurvRicErrStd'][0,:], label='R')[0]
    ax.plot(mcStats['time'], 3*mcStats['relVelCurvRicErrStd'][1,:], label='I')[0]
    ax.plot(mcStats['time'], 3*mcStats['relVelCurvRicErrStd'][2,:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Velocity Error')
    plt.grid()
       
    fig_curvRic_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    fig_curvRic_plt = plt.figure()
    plt.suptitle(r"Curvilinear RIC Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], mcStats['relPosCurvRicRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    plt.title('Relative Curvilinear RIC Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], mcStats['relVelCurvRicRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    plt.title('Relative Curvilinear RIC Velocity RMSE')
    plt.grid()
       
    fig_curvRic_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicRMSE_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DoeStats_Plot(mcStats, path, tag):
    """
    Delta Keplerian Orbit element statistics

    """
    
    fig_doe_plt = plt.figure()
    plt.suptitle(r"Differential Keplerian Orbit Element 3$\sigma$ Errors")
    
    ax = plt.subplot(3,2,1)
    ax.plot(mcStats['time'], 3*mcStats['doeErrStd'][0,:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta a$ [km]")
    plt.title('Delta Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(mcStats['time'], 3*mcStats['doeErrStd'][1,:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta e$")
    plt.title('Delta Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['doeErrStd'][2,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta i$ [deg]")
    plt.title('Delta Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['doeErrStd'][3,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta RAAN$ [deg]")
    plt.title('Delta RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['doeErrStd'][4,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta\omega$ [deg]")
    plt.title('Delta Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['doeErrStd'][5,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta M$ [deg]")
    plt.title('Delta Mean Anomaly')
    plt.grid()
       
    fig_doe_plt.tight_layout()
    
    fullFigPath = path + r"\doeStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefInrStats_Plot(mcStats, path, tag):
    """
    Chief inertial position and velocity statistics

    """
    
    fig_inr_plt = plt.figure()
    plt.suptitle(r"Chief Inertial 3$\sigma$ Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], 3*mcStats['posRsoErrStd'][0,:], label='X')[0]
    ax.plot(mcStats['time'], 3*mcStats['posRsoErrStd'][1,:], label='Y')[0]
    ax.plot(mcStats['time'], 3*mcStats['posRsoErrStd'][2,:], label='Z')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Chief Inertial Position Error')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], 3*mcStats['velRsoErrStd'][0,:], label='X')[0]
    ax.plot(mcStats['time'], 3*mcStats['velRsoErrStd'][1,:], label='Y')[0]
    ax.plot(mcStats['time'], 3*mcStats['velRsoErrStd'][2,:], label='Z')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Chief Inertial Velocity Error')
    plt.grid()
       
    fig_inr_plt.tight_layout()
    
    fullFigPath = path + r"\chiefInrStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    fig_inr_plt = plt.figure()
    plt.suptitle(r"Chief Inertial Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], mcStats['posRsoRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    plt.title('Chief Inertial Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], mcStats['velRsoRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    plt.title('Chief Inertial Velocity RMSE')
    plt.grid()
       
    fig_inr_plt.tight_layout()
    
    fullFigPath = path + r"\chiefInrRMSE_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefOeStats_Plot(mcStats, path, tag):
    """
    Chief Keplerian Element statistics

    """
    
    fig_oe_plt = plt.figure()
    plt.suptitle(r"Chief Orbit Element 3$\sigma$ Errors")
    
    ax = plt.subplot(3,2,1)
    ax.plot(mcStats['time'], 3*mcStats['oeRsoErrStd'][0,:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    plt.title('Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(mcStats['time'], 3*mcStats['oeRsoErrStd'][1,:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    plt.title('Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeRsoErrStd'][2,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    plt.title('Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeRsoErrStd'][3,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("RAAN [deg]")
    plt.title('RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeRsoErrStd'][4,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    plt.title('Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeRsoErrStd'][5,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    plt.title('Mean Anomaly')
    plt.grid()
       
    fig_oe_plt.tight_layout()
    
    fullFigPath = path + r"\chiefOeStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyInrStats_Plot(mcStats, path, tag):
    """
    Deputy inertial position and velocity statistics

    """
    
    fig_inr_plt = plt.figure()
    plt.suptitle(r"Deputy Inertial 3$\sigma$ Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], 3*mcStats['posVehErrStd'][0,:], label='X')[0]
    ax.plot(mcStats['time'], 3*mcStats['posVehErrStd'][1,:], label='Y')[0]
    ax.plot(mcStats['time'], 3*mcStats['posVehErrStd'][2,:], label='Z')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Deputy Inertial Position Error')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], 3*mcStats['velVehErrStd'][0,:], label='X')[0]
    ax.plot(mcStats['time'], 3*mcStats['velVehErrStd'][1,:], label='Y')[0]
    ax.plot(mcStats['time'], 3*mcStats['velVehErrStd'][2,:], label='Z')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Deputy Inertial Velocity Error')
    plt.grid()
       
    fig_inr_plt.tight_layout()
    
    fullFigPath = path + r"\deputyInrStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
    fig_inr_plt = plt.figure()
    plt.suptitle(r"Deputy Inertial Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    ax.plot(mcStats['time'], mcStats['posVehRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    plt.title('Deputy Inertial Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    ax.plot(mcStats['time'], mcStats['velVehRMSErr'])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    plt.title('Deputy Inertial Velocity RMSE')
    plt.grid()
       
    fig_inr_plt.tight_layout()
    
    fullFigPath = path + r"\deputyInrRMSE_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyOeStats_Plot(mcStats, path, tag):
    """
    Deputy Keplerian Element statistics

    """
    
    fig_oe_plt = plt.figure()
    plt.suptitle(r"Deputy Orbit Element 3$\sigma$ Errors")
    
    ax = plt.subplot(3,2,1)
    ax.plot(mcStats['time'], 3*mcStats['oeVehErrStd'][0,:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    plt.title('Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(mcStats['time'], 3*mcStats['oeVehErrStd'][1,:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    plt.title('Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeVehErrStd'][2,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    plt.title('Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeVehErrStd'][3,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("RAAN [deg]")
    plt.title('RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeVehErrStd'][4,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    plt.title('Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(mcStats['time'], 3*np.rad2deg(mcStats['oeVehErrStd'][5,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    plt.title('Mean Anomaly')
    plt.grid()
       
    fig_oe_plt.tight_layout()
    
    fullFigPath = path + r"\deputyOeStd_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ConsistencyStats_Plot(mcStats, path, tag):
    """
    Filter Consistency Stats Plots

    """
    
    fig_consistency_plt = plt.figure()
    
    ax1 = plt.subplot(2,1,1)
    ax1.plot(mcStats['time'], mcStats['avgNrmInvnSqr'])[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("NIS")
    plt.title('Average Normalized Innovations Squared')
    plt.grid()
    
    ax2 = plt.subplot(2,1,2)
    ax2.plot(mcStats['time'], mcStats['avgRsoNrmEstErrSqr'], label = 'Chief')[0]
    ax2.plot(mcStats['time'], mcStats['avgVehNrmEstErrSqr'], label = 'Deputy')[0]
    ax2.plot(mcStats['time'], mcStats['avgRelNrmEstErrSqr'], label = 'Relative')[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("NEES")
    ax2.legend()
    plt.title('Average Normalized Estimation Error Squared')
    plt.grid()
       
    fig_consistency_plt.tight_layout()
    
    fullFigPath = path + r"\consistencyStats_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Rectilinear RIC position and velocity statistics comparison

    """
    
    fig_rectRic_plt = plt.figure()
    plt.suptitle(r"Rectilinear RIC Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['relPosRectRicRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['relVelRectRicRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Velocity RMSE')
    plt.grid()
       
    fig_rectRic_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicRMSEComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Curvilinear RIC position and velocity statistics comparison

    """
    
    fig_curvRic_plt = plt.figure()
    plt.suptitle(r"Curvilinear RIC Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['relPosCurvRicRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['relVelCurvRicRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Velocity RMSE')
    plt.grid()
       
    fig_curvRic_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicRMSEComp_" + tag + r".png"
    plt.savefig(fullFigPath)    
    
def DoeStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Delta Keplerian Orbit element statistics comparison

    """
    
    fig_doe_plt = plt.figure()
    plt.suptitle(r"Differential Keplerian Orbit Element 3$\sigma$ Errors")
    
    ax = plt.subplot(3,2,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*mcStats[ii]['doeErrStd'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta a$ [km]")
    ax.legend()
    plt.title('Delta Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*mcStats[ii]['doeErrStd'][1,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta e$")
    ax.legend()
    plt.title('Delta Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['doeErrStd'][2,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta i$ [deg]")
    ax.legend()
    plt.title('Delta Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['doeErrStd'][3,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta RAAN$ [deg]")
    ax.legend()
    plt.title('Delta RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['doeErrStd'][4,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta\omega$ [deg]")
    ax.legend()
    plt.title('Delta Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['doeErrStd'][5,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta M$ [deg]")
    ax.legend()
    plt.title('Delta Mean Anomaly')
    plt.grid()
       
    fig_doe_plt.tight_layout()
    
    fullFigPath = path + r"\doeStdComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefInrStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Chief inertial position and velocity statistics comparison

    """
    
    fig_inr_plt = plt.figure()
    plt.suptitle(r"Chief Inertial Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['posRsoRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Chief Inertial Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['posRsoRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Chief Inertial Velocity RMSE')
    plt.grid()
       
    fig_inr_plt.tight_layout()
    
    fullFigPath = path + r"\chiefInrRMSEComp_" + tag + r".png"
    plt.savefig(fullFigPath)
  
def ChiefOeStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Chief Keplerian Element statistics comparison

    """
    
    fig_oe_plt = plt.figure()
    plt.suptitle(r"Chief Orbit Element 3$\sigma$ Errors")
    
    ax = plt.subplot(3,2,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*mcStats[ii]['oeRsoErrStd'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    ax.legend()
    plt.title('Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*mcStats[ii]['oeRsoErrStd'][1,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    ax.legend()
    plt.title('Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeRsoErrStd'][2,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    ax.legend()
    plt.title('Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeRsoErrStd'][3,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("RAAN [deg]")
    ax.legend()
    plt.title('RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeRsoErrStd'][4,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    ax.legend()
    plt.title('Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeRsoErrStd'][5,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    ax.legend()
    plt.title('Mean Anomaly')
    plt.grid()
       
    fig_oe_plt.tight_layout()
    
    fullFigPath = path + r"\chiefOeStdComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyInrStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Deputy inertial position and velocity statistics comparison

    """
    
    fig_inr_plt = plt.figure()
    plt.suptitle(r"Deputy Inertial Root Mean Square Errors")
    
    ax = plt.subplot(2,1,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['posVehRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Deputy Inertial Position RMSE')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['posVehRMSErr'], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Deputy Inertial Velocity RMSE')
    plt.grid()
       
    fig_inr_plt.tight_layout()
    
    fullFigPath = path + r"\deputyInrRMSEComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyOeStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Deputy Keplerian Element statistics comparison

    """
    
    fig_oe_plt = plt.figure()
    plt.suptitle(r"Deputy Orbit Element 3$\sigma$ Errors")
    
    ax = plt.subplot(3,2,1)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*mcStats[ii]['oeVehErrStd'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    ax.legend()
    plt.title('Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*mcStats[ii]['oeVehErrStd'][1,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    ax.legend()
    plt.title('Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeVehErrStd'][2,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    ax.legend()
    plt.title('Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeVehErrStd'][3,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("RAAN [deg]")
    ax.legend()
    plt.title('RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeVehErrStd'][4,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    ax.legend()
    plt.title('Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    # ax.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], 3*np.rad2deg(mcStats[ii]['oeVehErrStd'][5,:]), label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    ax.legend()
    plt.title('Mean Anomaly')
    plt.grid()
       
    fig_oe_plt.tight_layout()
    
    fullFigPath = path + r"\deputyOeStdComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ConsistencyStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Filter Consistency Stats comparison

    """
    
    fig_consistency_plt = plt.figure()
    
    ax1 = plt.subplot(4,1,1)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax1.plot(mcStats[ii]['time'], mcStats[ii]['avgNrmInvnSqr'], label=mcStr[ii])[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("NIS")
    ax1.legend()
    plt.title('Average Normalized Innovations Squared')
    plt.grid()
    
    ax2 = plt.subplot(4,1,2)
    # ax2.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax2.plot(mcStats[ii]['time'], mcStats[ii]['avgRsoNrmEstErrSqr'], label=mcStr[ii])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Chief NEES")
    ax2.legend()
    plt.title('Average Normalized Estimation Error Squared')
    plt.grid()
    
    ax3 = plt.subplot(4,1,3)
    # ax3.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax3.plot(mcStats[ii]['time'], mcStats[ii]['avgVehNrmEstErrSqr'], label=mcStr[ii])[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Deputy NEES")
    ax3.legend()
    plt.grid()
    
    ax4 = plt.subplot(4,1,4)
    # ax4.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax4.plot(mcStats[ii]['time'], mcStats[ii]['avgRelNrmEstErrSqr'], label=mcStr[ii])[0]
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Relative NEES")
    ax4.legend()
    plt.grid()
       
    fig_consistency_plt.tight_layout()
    
    fullFigPath = path + r"\consistencyStatsComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    