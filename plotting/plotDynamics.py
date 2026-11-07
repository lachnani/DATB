# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation

from plotting import plotGeneric as pg 

def visAll(log, path, tag, settings):
    """
    Plot all available plots and animations

    """
    plotAll(log, path, tag, settings)
    animAll(log, path, tag, settings)
    
def plotAll(log, path, tag, settings):
    """
    Plot all available plots

    """
    Fixed3D_Plot(log, path, tag)
    RectRicProj_Plot(log, path, tag)
    RectRic_Plot(log, path, tag)
    AccelDv_Plot(log, path, tag)
    if (settings["formation"]["relStates"] == True):
        ChiefOe_Plot(log, path, tag)
        ChiefEe_Plot(log, path, tag)
        CurvRicProj_Plot(log, path, tag)
        CurvRic_Plot(log, path, tag)
        Doe_Plot(log, path, tag)
        Dee_Plot(log, path, tag)
        RoePhase_Plot(log, path, tag)
        RectClroe_Plot(log, path, tag)
        CurvClroe_Plot(log, path, tag)
    if (settings["formation"]["environments"] == True):
        Environment_Plot(log, path, tag)
        SunMoonPos_Plot(log, path, tag)
    if (settings["formation"]["measurements"] == True):
        MeasParam_Plot(log, path, tag)
    
def animAll(log, path, tag, settings):
    """
    Plot all available animations

    """
    Fixed3D_Animation(log, path, tag)
    RectRicProj_Animation(log, path, tag)
    if (settings["formation"]["relStates"] == True):
        CurvRicProj_Animation(log, path, tag)


def Fixed3D_Plot(log, path, tag):
    """
    3D Plot from a fixed vantage point

    """
    
    fig_fixed3d_plt = plt.figure()
    ax = fig_fixed3d_plt.add_subplot(projection="3d")
    origin = ax.scatter(0, 0, 0, color='b', linewidths = 3)
    truthLine = ax.plot(log.relPosRectRic[1,:], log.relPosRectRic[0,:], log.relPosRectRic[2,:], color='r')
    truthEnd = ax.scatter(log.relPosRectRic[1,log.i], log.relPosRectRic[0,log.i], log.relPosRectRic[2,log.i], color='r')
    # Common steps
    ax.scatter(0, 0, 0, color='b', linewidths = 3)
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel("In-Track [km]")
    ax.set_ylabel("Radial [km]")
    ax.set_zlabel("Cross-Track [km]")
    #plt.grid()
    plt.gca().invert_xaxis()
    ax.view_init(elev=30, azim=-40, roll=0)
    plt.title('Relative Trajectory in RIC')
    
    fullFigPath = path + r"\fixed3d_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def Fixed3D_Animation(log, path, tag):
    """
    3D Animation from a fixed vantage point

    """
    
    def update_lines(frame, origin, truthLine, truthEnd):
        # update the line plot:
        R = log.relPosRectRic[0,:frame]
        I = log.relPosRectRic[1,:frame]
        C = log.relPosRectRic[2,:frame]
        truthLine.set_data_3d(I, R, C)
        pltMin = np.min([0, R.min(), I.min(), C.min()])
        pltMax = np.max([0, R.max(), I.max(), C.max()])
        pltLim = np.max(np.abs([pltMin, pltMax]))
        ax.set_xlim(-1*pltLim, pltLim)
        ax.set_ylim(-1*pltLim, pltLim)
        ax.set_zlim(-1*pltLim, pltLim)
        plt.gca().invert_xaxis()
        truthEnd.set_offsets(np.array((I[-1], R[-1])))
        truthEnd.set_3d_properties(C[-1],'z')
        plt.title('Relative Trajectory in RIC')
        return (truthLine, truthEnd)
    
    
    fig_fixed3d_ani = plt.figure()
    ax = fig_fixed3d_ani.add_subplot(projection="3d")
    origin = ax.scatter(0, 0, 0, color='b', linewidths = 3)
    truthLine = ax.plot(log.relPosRectRic[1,0], log.relPosRectRic[0,0], log.relPosRectRic[2,0], color='r')[0]
    truthEnd = ax.scatter(log.relPosRectRic[1,0], log.relPosRectRic[0,0], log.relPosRectRic[2,0], color='r')
    # Common steps
    ax.scatter(0, 0, 0, color='b', linewidths = 3)
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_xlabel("In-Track [km]")
    ax.set_ylabel("Radial [km]")
    ax.set_zlabel("Cross-Track [km]")
    #plt.grid()
    plt.gca().invert_xaxis()
    ax.view_init(elev=30, azim=-40, roll=0)
    
    ani_fixed3d_ani = animation.FuncAnimation(
        fig_fixed3d_ani, update_lines, log.i, fargs=(origin, truthLine, truthEnd), interval=1)
    
    fullFigPath = path + r"\fixed3dAni_" + tag + r".mp4"
    ani_fixed3d_ani.save(fullFigPath, fps=30)
    
    
def RectRicProj_Plot(log, path, tag):
    """
    2D Rectilinear RIC plots

    """
    
    fig_rectRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    riLine = ax1.plot(log.relPosRectRic[1,:], log.relPosRectRic[0,:], color='r')[0]
    riEnd = ax1.scatter(log.relPosRectRic[1,log.i], log.relPosRectRic[0,log.i], color='r')
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosRectRic[1,:], log.relPosRectRic[2,:], color='r')[0]
    icEnd = ax2.scatter(log.relPosRectRic[1,log.i], log.relPosRectRic[2,log.i], color='r')
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosRectRic[0,:], log.relPosRectRic[2,:], color='r')[0]
    rcEnd = ax3.scatter(log.relPosRectRic[0,log.i], log.relPosRectRic[2,log.i], color='r')
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Cross-Track [km]")
    ax3.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    fig_rectRicProj_plt.tight_layout()
    
    fullFigPath = path + r"\rectRicProj_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectRicProj_Animation(log, path, tag):
    """
    2D Rectilinear RIC animation

    """
    
    def update_lines(frame, origin, riLine, riEnd, icLine, icEnd, rcLine, rcEnd):
        # update the line plot:
        R = log.relPosRectRic[0,:frame]
        I = log.relPosRectRic[1,:frame]
        C = log.relPosRectRic[2,:frame]
        riLine.set_data(I, R)
        icLine.set_data(I, C)
        rcLine.set_data(R, C)
        riEnd.set_offsets(np.array((I[-1], R[-1])))
        icEnd.set_offsets(np.array((I[-1], C[-1])))
        rcEnd.set_offsets(np.array((R[-1], C[-1])))
        ax1.relim()
        ax1.autoscale()
        ax1.set_aspect('equal', adjustable='datalim')
        ax2.relim()
        ax2.autoscale()
        ax2.set_aspect('equal', adjustable='datalim')
        ax3.relim()
        ax3.autoscale()
        ax3.set_aspect('equal', adjustable='datalim')
        plt.tight_layout()
        plt.suptitle('Relative Trajectory in RIC')
        return (riLine, riEnd, icLine, icEnd, rcLine, rcEnd)
    
    fig_rectRicProj_ani = plt.figure()
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    riLine = ax1.plot(log.relPosRectRic[1,0], log.relPosRectRic[0,0], color='r')[0]
    riEnd = ax1.scatter(log.relPosRectRic[1,0], log.relPosRectRic[0,0], color='r')
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosRectRic[1,0], log.relPosRectRic[2,0], color='r')[0]
    icEnd = ax2.scatter(log.relPosRectRic[1,0], log.relPosRectRic[2,0], color='r')
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosRectRic[0,0], log.relPosRectRic[2,0], color='r')[0]
    rcEnd = ax3.scatter(log.relPosRectRic[0,0], log.relPosRectRic[2,0], color='r')
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Cross-Track [km]")
    plt.grid()
    
    fig_rectRicProj_ani.tight_layout()
    
    fig_rectRicProj_ani = animation.FuncAnimation(
        fig_rectRicProj_ani, update_lines, log.i, fargs=(origin, riLine, riEnd, icLine, icEnd, rcLine, rcEnd), interval=1)
    
    fullFigPath = path + r"\rectRicProjAni_" + tag + r".mp4"
    fig_rectRicProj_ani.save(fullFigPath, fps=30)
    
def CurvRicProj_Plot(log, path, tag):
    """
    2D Curvilinear RIC plots

    """
    
    fig_curvRicProj_plt = plt.figure()
    plt.suptitle('Relative Trajectory in RIC')
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    riLine = ax1.plot(log.relPosCurvRic[1,:], log.relPosCurvRic[0,:], color='r')[0]
    riEnd = ax1.scatter(log.relPosCurvRic[1,log.i], log.relPosCurvRic[0,log.i], color='r')
    ax1.set_xlabel("Curvilinear In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosCurvRic[1,:], log.relPosCurvRic[2,:], color='r')[0]
    icEnd = ax2.scatter(log.relPosCurvRic[1,log.i], log.relPosCurvRic[2,log.i], color='r')
    ax2.set_xlabel("Curvilinear In-Track [km]")
    ax2.set_ylabel("Curvilinear Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosCurvRic[0,:], log.relPosCurvRic[2,:], color='r')[0]
    rcEnd = ax3.scatter(log.relPosCurvRic[0,log.i], log.relPosCurvRic[2,log.i], color='r')
    ax3.set_xlabel("Radial [km]")
    ax3.set_ylabel("Curvilinear Cross-Track [km]")
    ax3.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    fig_curvRicProj_plt.tight_layout()
    
    fullFigPath = path + r"\curvRicProj_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRicProj_Animation(log, path, tag):
    """
    2D Curvilinear RIC animation

    """
    
    def update_lines(frame, origin, riLine, riEnd, icLine, icEnd, rcLine, rcEnd):
        # update the line plot:
        R = log.relPosCurvRic[0,:frame]
        I = log.relPosCurvRic[1,:frame]
        C = log.relPosCurvRic[2,:frame]
        riLine.set_data(I, R)
        icLine.set_data(I, C)
        rcLine.set_data(R, C)
        riEnd.set_offsets(np.array((I[-1], R[-1])))
        icEnd.set_offsets(np.array((I[-1], C[-1])))
        rcEnd.set_offsets(np.array((R[-1], C[-1])))
        ax1.relim()
        ax1.autoscale()
        ax1.set_aspect('equal', adjustable='datalim')
        ax2.relim()
        ax2.autoscale()
        ax2.set_aspect('equal', adjustable='datalim')
        ax3.relim()
        ax3.autoscale()
        ax3.set_aspect('equal', adjustable='datalim')
        plt.tight_layout()
        plt.suptitle('Relative Trajectory in RIC')
        return (riLine, riEnd, icLine, icEnd, rcLine, rcEnd)
    
    fig_curvRicProj_ani = plt.figure()
    
    ax1 = plt.subplot(212)
    origin = ax1.scatter(0, 0, color='b', linewidths = 3)
    riLine = ax1.plot(log.relPosCurvRic[1,0], log.relPosCurvRic[0,0], color='r')[0]
    riEnd = ax1.scatter(log.relPosCurvRic[1,0], log.relPosCurvRic[0,0], color='r')
    ax1.set_xlabel("Curvilinear In-Track [km]")
    ax1.set_ylabel("Curvilinear Radial [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosCurvRic[1,0], log.relPosCurvRic[2,0], color='r')[0]
    icEnd = ax2.scatter(log.relPosCurvRic[1,0], log.relPosCurvRic[2,0], color='r')
    ax2.set_xlabel("Curvilinear In-Track [km]")
    ax2.set_ylabel("Curvilinear Cross-Track [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosCurvRic[0,0], log.relPosCurvRic[2,0], color='r')[0]
    rcEnd = ax3.scatter(log.relPosCurvRic[0,0], log.relPosCurvRic[2,0], color='r')
    ax3.set_xlabel("Curvilinear Radial [km]")
    ax3.set_ylabel("Curvilinear Cross-Track [km]")
    plt.grid()
    
    fig_curvRicProj_ani.tight_layout()
    
    fig_curvRicProj_ani = animation.FuncAnimation(
        fig_curvRicProj_ani, update_lines, log.i, fargs=(origin, riLine, riEnd, icLine, icEnd, rcLine, rcEnd), interval=1)
    
    fullFigPath = path + r"\curvRicProjAni_" + tag + r".mp4"
    fig_curvRicProj_ani.save(fullFigPath, fps=30)
    
def RectRic_Plot(log, path, tag):
    """
    Rectilinear RIC position and velocities

    """   
    pg.plotStack(
        log.simTime, 
        [log.relPosRectRic,log.relVelRectRic], 
        [pg.line_ric,pg.line_ric], 
        pg.ax_posVel, 
        ['Relative Rectilinear RIC Position','Relative Rectilinear RIC Velocity'], 
        path + r"\rectRic_" + tag + r".png")
    
def CurvRic_Plot(log, path, tag):
    """
    Curvilinear RIC position and velocities

    """
    pg.plotStack(
        log.simTime, 
        [log.relPosRectRic,log.relVelRectRic], 
        [pg.line_ric,pg.line_ric], 
        pg.ax_posVel, 
        ['Relative Curvilinear RIC Position','Relative Curvilinear RIC Velocity'], 
        path + r"\curvRic_" + tag + r".png")
          
def AccelDv_Plot(log, path, tag):
    """
    Acceleration and accumulated DV plots

    """
    
    km2m = 1000.0
    
    fig = plt.figure()
    
    ax1 = plt.subplot(3,1,1)
    R = ax1.plot(log.simTime, log.accThrRsoRic[0,:]*km2m, label='R')[0]
    I = ax1.plot(log.simTime, log.accThrRsoRic[1,:]*km2m, label='I')[0]
    C = ax1.plot(log.simTime, log.accThrRsoRic[2,:]*km2m, label='C')[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Acceleration [m/$s^2$]")
    ax1.legend()
    plt.title('RSO Acceleration in RIC')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    R = ax2.plot(log.simTime, log.accThrVehRic[0,:]*km2m, label='R')[0]
    I = ax2.plot(log.simTime, log.accThrVehRic[1,:]*km2m, label='I')[0]
    C = ax2.plot(log.simTime, log.accThrVehRic[2,:]*km2m, label='C')[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Acceleration [m/$s^2$]")
    ax2.legend()
    plt.title('Vehicle Acceleration in RIC')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    rso = ax3.plot(log.simTime, log.dvTotRso*km2m, color='b', label='RSO')[0]
    veh = ax3.plot(log.simTime, log.dvTotVeh*km2m, color='r', label='Veh')[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("DV [m/s]")
    ax3.legend()
    plt.title('Accumulated Delta-V')
    plt.grid()
       
    fig.tight_layout()
    
    fullFigPath = path + r"\accelDv_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefOe_Plot(log, path, tag):
    """
    Chief Keplerian Elements plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.oeRso, 
        pg.ax_oe, 
        pg.title_oe, 
        path + r"\chiefOe_" + tag + r".png",
        suptitle='Chief Keplerian Orbit Elements')
    
def ChiefEe_Plot(log, path, tag):
    """
    Chief Equinoctial Elements plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.eeRso, 
        pg.ax_ee, 
        pg.title_ee, 
        path + r"\chiefEe_" + tag + r".png",
        suptitle='Chief Equinoctial Elements')
    
def Doe_Plot(log, path, tag):
    """
    Delta Keplerian Orbit elements plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.doe, 
        pg.ax_doe, 
        pg.title_doe, 
        path + r"\doe_" + tag + r".png",
        suptitle='Differential Keplerian Orbit Elements')
    
def Dee_Plot(log, path, tag):
    """
    Differential Equinoctial Elements plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.dee, 
        pg.ax_dee, 
        pg.title_dee, 
        path + r"\dee_" + tag + r".png",
        suptitle='Differential Equinoctial Elements')
    
def RectClroe_Plot(log, path, tag):
    """
    Rectilinear CLROE plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.rectClroe, 
        pg.ax_clroe, 
        pg.title_clroe, 
        path + r"\rectClroe_" + tag + r".png",
        suptitle='Rectilinear CLROEs')
    
def CurvClroe_Plot(log, path, tag):
    """
    Curvilinear CLROE plot

    """
    pg.plotThreeByTwo(
        log.simTime, 
        log.curvClroe, 
        pg.ax_clroe, 
        pg.title_clroe, 
        path + r"\rectClroe_" + tag + r".png",
        suptitle='Curvilinear CLROEs')
    
def Environment_Plot(log, path, tag):
    """
    Environment plots

    """
    
    fig_environment_plt = plt.figure()
    
    ax1 = plt.subplot(3,1,1)
    earth = ax1.plot(log.simTime, np.rad2deg(log.losAngEarth), label='Earth')[0]
    sun = ax1.plot(log.simTime, np.rad2deg(log.losAngSun), label='Sun')[0]
    moon = ax1.plot(log.simTime, np.rad2deg(log.losAngMoon), label='Moon')[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("LOS Angles [deg]")
    ax1.legend()
    plt.title('Line of Sight Angles')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    vis = ax2.plot(log.simTime, np.rad2deg(log.visibility))[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Visibility Angle [deg]")
    plt.title('Chief Visibility')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    rso = ax3.plot(log.simTime, log.eclipseRso, color='b', label='RSO')[0]
    veh = ax3.plot(log.simTime, log.eclipseVeh, color='r', label='Veh')[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Eclipse Flag")
    ax3.legend()
    plt.title('Eclipse State')
    plt.grid()
       
    fig_environment_plt.tight_layout()
    
    fullFigPath = path + r"\environment_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def SunMoonPos_Plot(log, path, tag):
    """
    Sun and Moon Position Plots

    """
    pg.plotStack(
        log.simTime, 
        [log.sunPosEci,log.moonPosEci], 
        [pg.line_vec,pg.line_vec], 
        pg.ax_posPos, 
        ['Sun Position in ECI','Moon Position in ECI'], 
        path + r"\sunMoonPos_" + tag + r".png")
    
def MeasParam_Plot(log, path, tag):
    """
    Measurement Parameter Plots

    """
    
    fig_measParam_plt = plt.figure()
    
    ax1 = plt.subplot(3,1,1)
    ax1.plot(log.simTime, log.measParams[0,:])[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Range [km]")
    plt.title('Relative Range')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    ax2.plot(log.simTime, log.measParams[1,:])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Range Rate [km/s]")
    plt.title('Relative Range Rate')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    az = ax3.plot(log.simTime, np.rad2deg(log.measParams[2,:]), label='Az')[0]
    el = ax3.plot(log.simTime, np.rad2deg(log.measParams[3,:]), label='El')[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("LOS Angles [deg]")
    ax3.legend()
    plt.title('Line of Sight Angles')
    plt.grid()
       
    fig_measParam_plt.tight_layout()
    
    fullFigPath = path + r"\measParam_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RoePhase_Plot(log, path, tag):
    """
    D'Amico ROE Phase plots

    """
    
    fig_roePhase_plt = plt.figure()
    plt.suptitle("ROE Phase")
    
    ax1 = plt.subplot(1,2,1)
    de = ax1.plot(log.roe[2,:], log.roe[3,:], color='r')[0]
    deEnd = ax1.scatter(log.roe[2,log.i], log.roe[3,log.i], color='r')
    ax1.set_xlabel(r"$\delta$$e_x$")
    ax1.set_ylabel(r"$\delta$$e_y$")
    # ax1.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    ax2 = plt.subplot(1,2,2)
    de = ax2.plot(log.roe[4,:], log.roe[5,:], color='r')[0]
    deEnd = ax2.scatter(log.roe[4,log.i], log.roe[5,log.i], color='r')
    ax2.set_xlabel(r"$\delta$$i_x$")
    ax2.set_ylabel(r"$\delta$$i_y$")
    # ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    fig_roePhase_plt.tight_layout()
    
    fullFigPath = path + r"\roePhase_" + tag + r".png"
    plt.savefig(fullFigPath)