# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

from plotting import plotGeneric as pg

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
    RelCovStats_Plot(mcStats, path, tag)
    ChiefInrStats_Plot(mcStats, path, tag)
    ChiefOeStats_Plot(mcStats, path, tag)
    ChiefCovStats_Plot(mcStats, path, tag)
    DeputyInrStats_Plot(mcStats, path, tag)
    DeputyOeStats_Plot(mcStats, path, tag)
    DeputyCovStats_Plot(mcStats, path, tag)
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
    CovarianceStatsComp_Plot(mcStats, mcStr, path, tag)
    
def RectRicProj_Plot(mcLog, path, tag):
    """
    2D Rectilinear RIC plots

    """
    
    log = mcLog[0]
    
    fig_rectRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    ax1.plot(log.relPosRectRic[1,:], log.relPosRectRic[0,:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosRectRic[1,log.i], log.relPosRectRic[0,log.i], color='r')
    for i in range(len(mcLog)):
        ax1.plot(mcLog[i].relPosRectRic[1,:], mcLog[i].fswNavRelPosRectRic[0,:], lw=0.5, color='k')[0]
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
    for i in range(len(mcLog)):
        ax2.plot(mcLog[i].relPosRectRic[1,:], mcLog[i].fswNavRelPosRectRic[2,:], lw=0.5, color='k')[0]
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
    for i in range(len(mcLog)):
        ax3.plot(mcLog[i].relPosRectRic[0,:], mcLog[i].fswNavRelPosRectRic[2,:], lw=0.5, color='k')[0]
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
    ax1.plot(log.relPosCurvRic[1,:], log.relPosCurvRic[0,:], color='r', label='Truth')[0]
    ax1.scatter(log.relPosCurvRic[1,log.i], log.relPosCurvRic[0,log.i], color='r')
    for i in range(len(mcLog)):
        ax1.plot(mcLog[i].relPosCurvRic[1,:], mcLog[i].fswNavRelPosCurvRic[0,:], lw=0.5, color='k')[0]
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    ax2.plot(log.relPosCurvRic[1,:], log.relPosCurvRic[2,:], color='r', label='Truth')[0]
    ax2.scatter(log.relPosCurvRic[1,log.i], log.relPosCurvRic[2,log.i], color='r')
    for i in range(len(mcLog)):
        ax2.plot(mcLog[i].relPosCurvRic[1,:], mcLog[i].fswNavRelPosCurvRic[2,:], lw=0.5, color='k')[0]
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend()
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    ax3.plot(log.relPosCurvRic[0,:], log.relPosCurvRic[2,:], color='r', label='Truth')[0]
    ax3.scatter(log.relPosCurvRic[0,log.i], log.relPosCurvRic[2,log.i], color='r')
    for i in range(len(mcLog)):
        ax3.plot(mcLog[i].relPosCurvRic[0,:], mcLog[i].fswNavRelPosCurvRic[2,:], lw=0.5, color='k')[0]
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
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelPosRectRicErr[0,:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelPosRectRicErr[1,:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelPosRectRicErr[2,:], lw=0.5, color='C2')[0]
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
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelVelRectRicErr[0,:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelVelRectRicErr[1,:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelVelRectRicErr[2,:], lw=0.5, color='C2')[0]
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
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelPosCurvRicErr[0,:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelPosCurvRicErr[1,:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelPosCurvRicErr[2,:], lw=0.5, color='C2')[0]
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
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelVelCurvRicErr[0,:], lw=0.5, color='C0')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Radial [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,2)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelVelCurvRicErr[1,:], lw=0.5, color='C1')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("In-Track [km/s]")
    plt.grid()
    
    ax = plt.subplot(3,1,3)
    for i in range(len(mcLog)):
        ax.plot(mcLog[i].simTime, mcLog[i].fswNavRelVelCurvRicErr[2,:], lw=0.5, color='C2')[0]
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
    pg.plotStack(
        mcStats['time'], 
        [3*mcStats['relPosRectRicErrStd'],3*mcStats['relVelRectRicErrStd']], 
        [pg.line_ric,pg.line_ric], 
        pg.ax_posVel, 
        ['Relative Rectilinear RIC Position Error','Relative Rectilinear RIC Velocity Error'], 
        path + r"\rectRicStd_" + tag + r".png",
        suptitle=r"Rectilinear RIC 3$\sigma$ Errors")
    
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
    pg.plotStack(
        mcStats['time'], 
        [3*mcStats['relPosCurvRicErrStd'],3*mcStats['relVelCurvRicErrStd']], 
        [pg.line_ric,pg.line_ric], 
        pg.ax_posVel, 
        ['Relative Curvilinear RIC Position Error','Relative Curvilinear RIC Velocity Error'], 
        path + r"\curvRicStd_" + tag + r".png",
        suptitle=r"Curvilinear RIC 3$\sigma$ Errors")
    
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
    pg.plotThreeByTwo(
        mcStats['time'], 
        3*mcStats['doeErrStd'], 
        pg.ax_doe, 
        pg.title_doe, 
        path + r"\doeStd_" + tag + r".png",
        suptitle=r"Differential Keplerian Orbit Element 3$\sigma$ Errors")
    
def RelCovStats_Plot(mcStats, path, tag):
    """
    Relative RIC covariance statistics
    
    """
    
    fig = plt.figure()
    plt.suptitle('Relative RIC Covariance Statistics')
    
    ax = plt.subplot(2,2,1)
    ax.plot(mcStats['time'], np.sqrt(mcStats['avgRelCovDet'][0,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^3]")
    plt.title('Average Position Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(2,2,2)
    ax.plot(mcStats['time'], mcStats['avgRelCovTrace'][0,:], label='Trace')[0]
    ax.plot(mcStats['time'], mcStats['avgRelCovEig'][0,:], label='Max Eig')[0]
    ax.plot(mcStats['time'], mcStats['avgRelCovFro'][0,:], label='Frobenius Norm')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^2]")
    ax.legend()
    plt.title('Average Position Covariance Size')
    plt.grid()
    
    ax = plt.subplot(2,2,3)
    ax.plot(mcStats['time'], np.sqrt(mcStats['avgRelCovDet'][1,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[(km/s)^3]")
    plt.title('Average Velocity Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(2,2,4)
    ax.plot(mcStats['time'], mcStats['avgRelCovTrace'][1,:], label='Trace')[0]
    ax.plot(mcStats['time'], mcStats['avgRelCovEig'][3,:], label='Max Eig')[0]
    ax.plot(mcStats['time'], mcStats['avgRelCovFro'][1,:], label='Frobenius Norm')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[(km/s)^2]")
    ax.legend()
    plt.title('Average Velocity Covariance Size')
    plt.grid()
       
    fig.tight_layout()
    
    fullFigPath = path + r"\relCovStats_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefInrStats_Plot(mcStats, path, tag):
    """
    Chief inertial position and velocity statistics

    """
    pg.plotStack(
        mcStats['time'], 
        [3*mcStats['posRsoErrStd'],3*mcStats['velRsoErrStd']], 
        [pg.line_vec,pg.line_vec], 
        pg.ax_posVel, 
        ['Chief Inertial Position Error','Chief Inertial Velocity Error'], 
        path + r"\chiefInrStd_" + tag + r".png",
        suptitle=r"Chief Inertial 3$\sigma$ Errors")
    
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
    pg.plotThreeByTwo(
        mcStats['time'], 
        3*mcStats['oeRsoErrStd'], 
        pg.ax_oe, 
        pg.title_oe, 
        path + r"\chiefOeStd_" + tag + r".png",
        suptitle=r"Chief Orbit Element 3$\sigma$ Errors")
    
def ChiefCovStats_Plot(mcStats, path, tag):
    """
    Chief covariance statistics
    
    """
    
    fig = plt.figure()
    plt.suptitle('Chief Covariance Statistics')
    
    ax = plt.subplot(2,2,1)
    ax.plot(mcStats['time'], np.sqrt(mcStats['avgRsoCovDet'][0,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^3]")
    plt.title('Average Position Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(2,2,2)
    ax.plot(mcStats['time'], mcStats['avgRsoCovTrace'][0,:], label='Trace')[0]
    ax.plot(mcStats['time'], mcStats['avgRsoCovEig'][0,:], label='Max Eig')[0]
    ax.plot(mcStats['time'], mcStats['avgRsoCovFro'][0,:], label='Frobenius Norm')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^2]")
    ax.legend()
    plt.title('Average Position Covariance Size')
    plt.grid()
    
    ax = plt.subplot(2,2,3)
    ax.plot(mcStats['time'], np.sqrt(mcStats['avgRsoCovDet'][1,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[(km/s)^3]")
    plt.title('Average Velocity Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(2,2,4)
    ax.plot(mcStats['time'], mcStats['avgRsoCovTrace'][1,:], label='Trace')[0]
    ax.plot(mcStats['time'], mcStats['avgRsoCovEig'][3,:], label='Max Eig')[0]
    ax.plot(mcStats['time'], mcStats['avgRsoCovFro'][1,:], label='Frobenius Norm')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[(km/s)^2]")
    ax.legend()
    plt.title('Average Velocity Covariance Size')
    plt.grid()
       
    fig.tight_layout()
    
    fullFigPath = path + r"\chiefCovStats_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def DeputyInrStats_Plot(mcStats, path, tag):
    """
    Deputy inertial position and velocity statistics

    """
    pg.plotStack(
        mcStats['time'], 
        [3*mcStats['posVehErrStd'],3*mcStats['velVehErrStd']], 
        [pg.line_vec,pg.line_vec], 
        pg.ax_posVel, 
        ['Deputy Inertial Position Error','Deputy Inertial Velocity Error'], 
        path + r"\deputyInrStd_" + tag + r".png",
        suptitle=r"Deputy Inertial 3$\sigma$ Errors")
    
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
    pg.plotThreeByTwo(
        mcStats['time'], 
        3*mcStats['oeVehErrStd'], 
        pg.ax_oe, 
        pg.title_oe, 
        path + r"\deputyOeStd_" + tag + r".png",
        suptitle=r"Deputy Orbit Element 3$\sigma$ Errors")
    
def DeputyCovStats_Plot(mcStats, path, tag):
    """
    Deputy covariance statistics
    
    """
    
    fig = plt.figure()
    plt.suptitle('Deputy Covariance Statistics')
    
    ax = plt.subplot(2,2,1)
    ax.plot(mcStats['time'], np.sqrt(mcStats['avgVehCovDet'][0,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^3]")
    plt.title('Average Position Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(2,2,2)
    ax.plot(mcStats['time'], mcStats['avgVehCovTrace'][0,:], label='Trace')[0]
    ax.plot(mcStats['time'], mcStats['avgVehCovEig'][0,:], label='Max Eig')[0]
    ax.plot(mcStats['time'], mcStats['avgVehCovFro'][0,:], label='Frobenius Norm')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^2]")
    ax.legend()
    plt.title('Average Position Covariance Size')
    plt.grid()
    
    ax = plt.subplot(2,2,3)
    ax.plot(mcStats['time'], np.sqrt(mcStats['avgVehCovDet'][1,:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[(km/s)^3]")
    plt.title('Average Velocity Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(2,2,4)
    ax.plot(mcStats['time'], mcStats['avgVehCovTrace'][1,:], label='Trace')[0]
    ax.plot(mcStats['time'], mcStats['avgVehCovEig'][3,:], label='Max Eig')[0]
    ax.plot(mcStats['time'], mcStats['avgVehCovFro'][1,:], label='Frobenius Norm')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[(km/s)^2]")
    ax.legend()
    plt.title('Average Velocity Covariance Size')
    plt.grid()
       
    fig.tight_layout()
    
    fullFigPath = path + r"\deputyCovStats_" + tag + r".png"
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
        ax.plot(mcStats[ii]['time'], mcStats[ii]['velRsoRMSErr'], label=mcStr[ii])[0]
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
        ax.plot(mcStats[ii]['time'], mcStats[ii]['velVehRMSErr'], label=mcStr[ii])[0]
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
    
def CovarianceStatsComp_Plot(mcStats, mcStr, path, tag):
    """
    Covariance Stats comparison

    """
    
    fig = plt.figure()
    
    ax = plt.subplot(3,2,1)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['avgRsoCovDet'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^3]")
    ax.legend()
    plt.title('Average Chief Position Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['avgRsoCovTrace'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^2]")
    ax.legend()
    plt.title('Average Chief Position Trace')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['avgVehCovDet'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^3]")
    ax.legend()
    plt.title('Average Deputy Position Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['avgVehCovTrace'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^2]")
    ax.legend()
    plt.title('Average Deputy Position Trace')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['avgRelCovDet'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^3]")
    ax.legend()
    plt.title('Average Relative RIC Position Sqrt Determinant')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    # ax1.set_prop_cycle(linestyle=["-", "--", "-.", ":"])
    for ii in range(len(mcStats)):
        ax.plot(mcStats[ii]['time'], mcStats[ii]['avgRelCovTrace'][0,:], label=mcStr[ii])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("[km^2]")
    ax.legend()
    plt.title('Average Relative RIC Position Trace')
    plt.grid()
       
    fig.tight_layout()
    
    fullFigPath = path + r"\covarianceStatsComp_" + tag + r".png"
    plt.savefig(fullFigPath)
    