# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation

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
    truthLine = ax.plot(log.relPosRectRic.loc['I',:], log.relPosRectRic.loc['R',:], log.relPosRectRic.loc['C',:], color='r')
    truthEnd = ax.scatter(log.relPosRectRic.loc['I',log.i], log.relPosRectRic.loc['R',log.i], log.relPosRectRic.loc['C',log.i], color='r')
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
        R = log.relPosRectRic.loc['R',:frame].to_numpy()
        I = log.relPosRectRic.loc['I',:frame].to_numpy()
        C = log.relPosRectRic.loc['C',:frame].to_numpy()
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
    truthLine = ax.plot(log.relPosRectRic.loc['I',0], log.relPosRectRic.loc['R',0], log.relPosRectRic.loc['C',0], color='r')[0]
    truthEnd = ax.scatter(log.relPosRectRic.loc['I',0], log.relPosRectRic.loc['R',0], log.relPosRectRic.loc['C',0], color='r')
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
    riLine = ax1.plot(log.relPosRectRic.loc['I',:], log.relPosRectRic.loc['R',:], color='r')[0]
    riEnd = ax1.scatter(log.relPosRectRic.loc['I',log.i], log.relPosRectRic.loc['R',log.i], color='r')
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosRectRic.loc['I',:], log.relPosRectRic.loc['C',:], color='r')[0]
    icEnd = ax2.scatter(log.relPosRectRic.loc['I',log.i], log.relPosRectRic.loc['C',log.i], color='r')
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosRectRic.loc['R',:], log.relPosRectRic.loc['C',:], color='r')[0]
    rcEnd = ax3.scatter(log.relPosRectRic.loc['R',log.i], log.relPosRectRic.loc['C',log.i], color='r')
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
        R = log.relPosRectRic.loc['R',:frame].to_numpy()
        I = log.relPosRectRic.loc['I',:frame].to_numpy()
        C = log.relPosRectRic.loc['C',:frame].to_numpy()
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
    riLine = ax1.plot(log.relPosRectRic.loc['I',0], log.relPosRectRic.loc['R',0], color='r')[0]
    riEnd = ax1.scatter(log.relPosRectRic.loc['I',0], log.relPosRectRic.loc['R',0], color='r')
    ax1.set_xlabel("In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosRectRic.loc['I',0], log.relPosRectRic.loc['C',0], color='r')[0]
    icEnd = ax2.scatter(log.relPosRectRic.loc['I',0], log.relPosRectRic.loc['C',0], color='r')
    ax2.set_xlabel("In-Track [km]")
    ax2.set_ylabel("Cross-Track [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosRectRic.loc['R',0], log.relPosRectRic.loc['C',0], color='r')[0]
    rcEnd = ax3.scatter(log.relPosRectRic.loc['R',0], log.relPosRectRic.loc['C',0], color='r')
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
    riLine = ax1.plot(log.relPosCurvRic.loc['I',:], log.relPosCurvRic.loc['R',:], color='r')[0]
    riEnd = ax1.scatter(log.relPosCurvRic.loc['I',log.i], log.relPosCurvRic.loc['R',log.i], color='r')
    ax1.set_xlabel("Curvilinear In-Track [km]")
    ax1.set_ylabel("Radial [km]")
    ax1.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosCurvRic.loc['I',:], log.relPosCurvRic.loc['C',:], color='r')[0]
    icEnd = ax2.scatter(log.relPosCurvRic.loc['I',log.i], log.relPosCurvRic.loc['C',log.i], color='r')
    ax2.set_xlabel("Curvilinear In-Track [km]")
    ax2.set_ylabel("Curvilinear Cross-Track [km]")
    ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosCurvRic.loc['R',:], log.relPosCurvRic.loc['C',:], color='r')[0]
    rcEnd = ax3.scatter(log.relPosCurvRic.loc['R',log.i], log.relPosCurvRic.loc['C',log.i], color='r')
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
        R = log.relPosCurvRic.loc['R',:frame].to_numpy()
        I = log.relPosCurvRic.loc['I',:frame].to_numpy()
        C = log.relPosCurvRic.loc['C',:frame].to_numpy()
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
    riLine = ax1.plot(log.relPosCurvRic.loc['I',0], log.relPosCurvRic.loc['R',0], color='r')[0]
    riEnd = ax1.scatter(log.relPosCurvRic.loc['I',0], log.relPosCurvRic.loc['R',0], color='r')
    ax1.set_xlabel("Curvilinear In-Track [km]")
    ax1.set_ylabel("Curvilinear Radial [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax2 = plt.subplot(221)
    origin = ax2.scatter(0, 0, color='b', linewidths = 3)
    icLine = ax2.plot(log.relPosCurvRic.loc['I',0], log.relPosCurvRic.loc['C',0], color='r')[0]
    icEnd = ax2.scatter(log.relPosCurvRic.loc['I',0], log.relPosCurvRic.loc['C',0], color='r')
    ax2.set_xlabel("Curvilinear In-Track [km]")
    ax2.set_ylabel("Curvilinear Cross-Track [km]")
    plt.grid()
    plt.gca().invert_xaxis()
    
    ax3 = plt.subplot(222)
    origin = ax3.scatter(0, 0, color='b', linewidths = 3)
    rcLine = ax3.plot(log.relPosCurvRic.loc['R',0], log.relPosCurvRic.loc['C',0], color='r')[0]
    rcEnd = ax3.scatter(log.relPosCurvRic.loc['R',0], log.relPosCurvRic.loc['C',0], color='r')
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
    
    fig_rectRic_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.relPosRectRic.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.relPosRectRic.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.relPosRectRic.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Position')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.relVelRectRic.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.relVelRectRic.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.relVelRectRic.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Rectilinear RIC Velocity')
    plt.grid()
       
    fig_rectRic_plt.tight_layout()
    
    fullFigPath = path + r"\rectRic_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvRic_Plot(log, path, tag):
    """
    Curvilinear RIC position and velocities

    """
    
    fig_curvRic_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.relPosCurvRic.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.relPosCurvRic.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.relPosCurvRic.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Position')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.relVelCurvRic.loc['R',:], label='R')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.relVelCurvRic.loc['I',:], label='I')[0]
    X = ax.plot(log.time.loc['simTime',:], log.relVelCurvRic.loc['C',:], label='C')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.legend()
    plt.title('Relative Curvilinear RIC Velocity')
    plt.grid()
       
    fig_curvRic_plt.tight_layout()
    
    fullFigPath = path + r"\curvRic_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def AccelDv_Plot(log, path, tag):
    """
    Acceleration and accumulated DV plots

    """
    
    km2m = 1000.0
    
    fig_accelDv_plt = plt.figure()
    
    ax1 = plt.subplot(3,1,1)
    R = ax1.plot(log.time.loc['simTime',:], log.accThrRsoRic.loc['R',:]*km2m, label='R')[0]
    I = ax1.plot(log.time.loc['simTime',:], log.accThrRsoRic.loc['I',:]*km2m, label='I')[0]
    C = ax1.plot(log.time.loc['simTime',:], log.accThrRsoRic.loc['C',:]*km2m, label='C')[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Acceleration [m/$s^2$]")
    ax1.legend()
    plt.title('RSO Acceleration in RIC')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    R = ax2.plot(log.time.loc['simTime',:], log.accThrVehRic.loc['R',:]*km2m, label='R')[0]
    I = ax2.plot(log.time.loc['simTime',:], log.accThrVehRic.loc['I',:]*km2m, label='I')[0]
    C = ax2.plot(log.time.loc['simTime',:], log.accThrVehRic.loc['C',:]*km2m, label='C')[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Acceleration [m/$s^2$]")
    ax2.legend()
    plt.title('Vehicle Acceleration in RIC')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    rso = ax3.plot(log.time.loc['simTime',:], log.dvTotRso.loc['DV',:]*km2m, color='b', label='RSO')[0]
    veh = ax3.plot(log.time.loc['simTime',:], log.dvTotVeh.loc['DV',:]*km2m, color='r', label='Veh')[0]
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("DV [m/s]")
    ax3.legend()
    plt.title('Accumulated Delta-V')
    plt.grid()
       
    fig_accelDv_plt.tight_layout()
    
    fullFigPath = path + r"\accelDv_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefOe_Plot(log, path, tag):
    """
    Chief Keplerian Elements plot

    """
    
    fig_chiefOe_plt = plt.figure()
    plt.suptitle('Chief Keplerian Orbit Elements')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.oeRso.loc['a',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    plt.title('Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], log.oeRso.loc['e',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("e")
    plt.title('Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['i',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("i [deg]")
    plt.title('Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['RAAN',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\Omega$ [deg]")
    plt.title('RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['argP',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\omega$ [deg]")
    plt.title('Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.oeRso.loc['M',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("M [deg]")
    plt.title('Mean Anomaly')
    plt.grid()
       
    fig_chiefOe_plt.tight_layout()
    
    fullFigPath = path + r"\chiefOe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def ChiefEe_Plot(log, path, tag):
    """
    Chief Equinoctial Elements plot

    """
    
    fig_chiefEe_plt = plt.figure()
    plt.suptitle('Chief Equinoctial Elements')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.eeRso.loc['a',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("a [km]")
    plt.title('Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.eeRso.loc['l',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("l [deg]")
    plt.title('Mean Longitude')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.eeRso.loc['P1',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("P1")
    plt.title('P1')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.eeRso.loc['P2',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("P2")
    plt.title('P2')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.eeRso.loc['Q1',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Q1")
    plt.title('Q1')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], log.eeRso.loc['Q2',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Q2")
    plt.title('Q2')
    plt.grid()
       
    fig_chiefEe_plt.tight_layout()
    
    fullFigPath = path + r"\chiefEe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def Doe_Plot(log, path, tag):
    """
    Delta Keplerian Orbit elements plot

    """
    
    fig_doe_plt = plt.figure()
    plt.suptitle('Differential Keplerian Orbit Elements')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.doe.loc['da',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta a$ [km]")
    plt.title('Delta Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], log.doe.loc['de',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta e$")
    plt.title('Delta Eccentricity')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.doe.loc['di',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta i$ [deg]")
    plt.title('Delta Inclination')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.doe.loc['dRAAN',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta RAAN$ [deg]")
    plt.title('Delta RAAN')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.doe.loc['dargP',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta\omega$ [deg]")
    plt.title('Delta Argument of Perigee')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.doe.loc['dM',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta M$ [deg]")
    plt.title('Delta Mean Anomaly')
    plt.grid()
       
    fig_doe_plt.tight_layout()
    
    fullFigPath = path + r"\doe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def Dee_Plot(log, path, tag):
    """
    Differential Equinoctial Elements plot

    """
    
    fig_dee_plt = plt.figure()
    plt.suptitle('Differential Equinoctial Elements')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.dee.loc['da',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta a [km]$")
    plt.title('Delta Semimajor-Axis')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.dee.loc['dl',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta l [deg]$")
    plt.title('Delta Mean Longitude')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.dee.loc['dP1',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta P1$")
    plt.title('Delta P1')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.dee.loc['dP2',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta P2$")
    plt.title('Delta P2')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.dee.loc['dQ1',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta Q1$")
    plt.title('Delta Q1')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], log.dee.loc['dQ2',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\delta Q2$")
    plt.title('Delta Q2')
    plt.grid()
       
    fig_dee_plt.tight_layout()
    
    fullFigPath = path + r"\dee_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def RectClroe_Plot(log, path, tag):
    """
    Rectilinear CLROE plot

    """
    
    fig_rectClroe_plt = plt.figure()
    plt.suptitle('Rectilinear CLROEs')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.rectClroe.loc['A0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$A_0$ [km]")
    plt.title('In-Plane Ellipse Size')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.rectClroe.loc['alpha',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\alpha$ [deg]")
    plt.title('In-Plane Ellipse Phase')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.rectClroe.loc['xOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$x_{off}$ [km]")
    plt.title('Radial Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.rectClroe.loc['yOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$y_{off}$ [km]")
    plt.title('In-Track Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.rectClroe.loc['B0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$B_0$ [km]")
    plt.title('Cross-Track Magnitude')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.rectClroe.loc['beta',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\beta$ [deg]")
    plt.title('Cross-Track Phase')
    plt.grid()
       
    fig_rectClroe_plt.tight_layout()
    
    fullFigPath = path + r"\rectClroe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def CurvClroe_Plot(log, path, tag):
    """
    Curvilinear CLROE plot

    """
    
    fig_curvClroe_plt = plt.figure()
    plt.suptitle('Curvilinear CLROEs')
    
    ax = plt.subplot(3,2,1)
    ax.plot(log.time.loc['simTime',:], log.curvClroe.loc['A0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$A_0$ [km]")
    plt.title('In-Plane Ellipse Size')
    plt.grid()
    
    ax = plt.subplot(3,2,2)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.curvClroe.loc['alpha',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\alpha$ [deg]")
    plt.title('In-Plane Ellipse Phase')
    plt.grid()
    
    ax = plt.subplot(3,2,3)
    ax.plot(log.time.loc['simTime',:], log.curvClroe.loc['xOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$x_{off}$ [km]")
    plt.title('Radial Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,4)
    ax.plot(log.time.loc['simTime',:], log.curvClroe.loc['yOff',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$y_{off}$ [km]")
    plt.title('In-Track Offset')
    plt.grid()
    
    ax = plt.subplot(3,2,5)
    ax.plot(log.time.loc['simTime',:], log.curvClroe.loc['B0',:])[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$B_0$ [km]")
    plt.title('Cross-Track Magnitude')
    plt.grid()
    
    ax = plt.subplot(3,2,6)
    ax.plot(log.time.loc['simTime',:], np.rad2deg(log.curvClroe.loc['beta',:]))[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel(r"$\beta$ [deg]")
    plt.title('Cross-Track Phase')
    plt.grid()
       
    fig_curvClroe_plt.tight_layout()
    
    fullFigPath = path + r"\curvClroe_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def Environment_Plot(log, path, tag):
    """
    Environment plots

    """
    
    fig_environment_plt = plt.figure()
    
    ax1 = plt.subplot(2,1,1)
    earth = ax1.plot(log.time.loc['simTime',:], np.rad2deg(log.losAngles.loc['Earth',:]), label='Earth')[0]
    sun = ax1.plot(log.time.loc['simTime',:], np.rad2deg(log.losAngles.loc['Sun',:]), label='Sun')[0]
    moon = ax1.plot(log.time.loc['simTime',:], np.rad2deg(log.losAngles.loc['Moon',:]), label='Moon')[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("LOS Angles [deg]")
    ax1.legend()
    plt.title('Line of Sight Angles')
    plt.grid()
    
    ax2 = plt.subplot(2,1,2)
    rso = ax2.plot(log.time.loc['simTime',:], log.eclipseRso.loc['Eclipse',:], color='b', label='RSO')[0]
    veh = ax2.plot(log.time.loc['simTime',:], log.eclipseVeh.loc['Eclipse',:], color='r', label='Veh')[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Eclipse Flag")
    ax2.legend()
    plt.title('Eclipse State')
    plt.grid()
       
    fig_environment_plt.tight_layout()
    
    fullFigPath = path + r"\environment_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def SunMoonPos_Plot(log, path, tag):
    """
    Sun and Moon Position Plots

    """
    
    fig_sunMoonPos_plt = plt.figure()
    
    ax = plt.subplot(2,1,1)
    X = ax.plot(log.time.loc['simTime',:], log.sunPosEci.loc['X',:], label='X')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.sunPosEci.loc['Y',:], label='Y')[0]
    Z = ax.plot(log.time.loc['simTime',:], log.sunPosEci.loc['Z',:], label='Z')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Sun Position in ECI')
    plt.grid()
    
    ax = plt.subplot(2,1,2)
    X = ax.plot(log.time.loc['simTime',:], log.moonPosEci.loc['X',:], label='X')[0]
    Y = ax.plot(log.time.loc['simTime',:], log.moonPosEci.loc['Y',:], label='Y')[0]
    Z = ax.plot(log.time.loc['simTime',:], log.moonPosEci.loc['Z',:], label='Z')[0]
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position [km]")
    ax.legend()
    plt.title('Moon Position in ECI')
    plt.grid()
       
    fig_sunMoonPos_plt.tight_layout()
    
    fullFigPath = path + r"\sunMoonPos_" + tag + r".png"
    plt.savefig(fullFigPath)
    
def MeasParam_Plot(log, path, tag):
    """
    Measurement Parameter Plots

    """
    
    fig_measParam_plt = plt.figure()
    
    ax1 = plt.subplot(3,1,1)
    ax1.plot(log.time.loc['simTime',:], log.measParams.loc['Rng',:])[0]
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Range [km]")
    plt.title('Relative Range')
    plt.grid()
    
    ax2 = plt.subplot(3,1,2)
    ax2.plot(log.time.loc['simTime',:], log.measParams.loc['Rng-Rate',:])[0]
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Range Rate [km/s]")
    plt.title('Relative Range Rate')
    plt.grid()
    
    ax3 = plt.subplot(3,1,3)
    az = ax3.plot(log.time.loc['simTime',:], np.rad2deg(log.measParams.loc['Az',:]), label='Az')[0]
    el = ax3.plot(log.time.loc['simTime',:], np.rad2deg(log.measParams.loc['El',:]), label='El')[0]
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
    plt.suptitle("D'Amico ROE Phase")
    
    ax1 = plt.subplot(1,3,1)
    aLambda = ax1.plot(log.roe.loc['dlambda',:], log.roe.loc['da',:], color='r')[0]
    aLambdaEnd = ax1.scatter(log.roe.loc['dlambda',log.i], log.roe.loc['da',log.i], color='r')
    ax1.set_xlabel(r"$\delta\lambda$")
    ax1.set_ylabel(r"$\delta$a")
    # ax1.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    ax2 = plt.subplot(1,3,2)
    de = ax2.plot(log.roe.loc['dex',:], log.roe.loc['dey',:], color='r')[0]
    deEnd = ax2.scatter(log.roe.loc['dex',log.i], log.roe.loc['dey',log.i], color='r')
    ax2.set_xlabel(r"$\delta$$e_x$")
    ax2.set_ylabel(r"$\delta$$e_y$")
    # ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    ax3 = plt.subplot(1,3,3)
    de = ax3.plot(log.roe.loc['dix',:], log.roe.loc['diy',:], color='r')[0]
    deEnd = ax3.scatter(log.roe.loc['dix',log.i], log.roe.loc['diy',log.i], color='r')
    ax3.set_xlabel(r"$\delta$$i_x$")
    ax3.set_ylabel(r"$\delta$$i_y$")
    # ax2.set_aspect('equal', adjustable='datalim')
    plt.grid()
    
    fig_roePhase_plt.tight_layout()
    
    fullFigPath = path + r"\roePhase_" + tag + r".png"
    plt.savefig(fullFigPath)