# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

def visAll(mcLog, mcStats, path, tag, settings):
    """
    Plot all available plots and animations

    """
    plotAll(mcLog, mcStats, path, tag, settings)
    
def plotAll(mcLog, mcStats, path, tag, settings):
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
    ChiefOeStats_Plot(mcStats, path, tag)
    DeputyOeStats_Plot(mcStats, path, tag)
    ConsistencyStats_Plot(mcStats, path, tag)
    

    
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
    
    fullFigPath = path + r"\rectRicStats_" + tag + r".png"
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
    
    fullFigPath = path + r"\rectRicStats_" + tag + r".png"
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
    
    fullFigPath = path + r"\chiefOeStats_" + tag + r".png"
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
    
    fullFigPath = path + r"\deputyOeStats_" + tag + r".png"
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
    ax2.plot(mcStats['time'], mcStats['avgNrmEstErrSqr'])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("NEES")
    plt.title('Average Normalized Estimation Error Squared')
    plt.grid()
       
    fig_consistency_plt.tight_layout()
    
    fullFigPath = path + r"\consistencyStats_" + tag + r".png"
    plt.savefig(fullFigPath)