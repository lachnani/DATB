# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:12:19 2026

@author: Hakim Lachnani
"""

import numpy as np

def computeMcStats(mcLog, numSims):
    print("MC: Accumulating Data ...")
    sim_len = len(mcLog[0].time.loc['simTime',:])
    mcRelPosRectRicErr = np.zeros((3,sim_len,numSims))
    mcRelVelRectRicErr = np.zeros((3,sim_len,numSims))
    mcRelPosCurvRicErr = np.zeros((3,sim_len,numSims))
    mcRelVelCurvRicErr = np.zeros((3,sim_len,numSims))
    mcDoeErr = np.zeros((6,sim_len,numSims))
    mcPosVehErr = np.zeros((3,sim_len,numSims))
    mcVelVehErr = np.zeros((3,sim_len,numSims))
    mcOeVehErr = np.zeros((6,sim_len,numSims))
    mcPosRsoErr = np.zeros((3,sim_len,numSims))
    mcVelRsoErr = np.zeros((3,sim_len,numSims))
    mcOeRsoErr = np.zeros((6,sim_len,numSims))
    mcNrmInvnSqr = np.zeros((sim_len,numSims))
    mcNrmEstErrSqr = np.zeros((sim_len,numSims))
    for simNum in range(numSims):
        mcRelPosRectRicErr[:,:,simNum] = mcLog[simNum].fswNavRelPosRectRicErr.to_numpy()
        mcRelVelRectRicErr[:,:,simNum] = mcLog[simNum].fswNavRelVelRectRicErr.to_numpy()
        mcRelPosCurvRicErr[:,:,simNum] = mcLog[simNum].fswNavRelPosCurvRicErr.to_numpy()
        mcRelVelCurvRicErr[:,:,simNum] = mcLog[simNum].fswNavRelVelCurvRicErr.to_numpy()
        mcDoeErr[:,:,simNum] = mcLog[simNum].fswNavDoeErr.to_numpy()
        mcPosVehErr[:,:,simNum] = mcLog[simNum].fswNavPosVehEciErr.to_numpy()
        mcVelVehErr[:,:,simNum] = mcLog[simNum].fswNavVelVehEciErr.to_numpy()
        mcOeVehErr[:,:,simNum] = mcLog[simNum].oeVeh.to_numpy() - mcLog[simNum].fswNavOeVeh.to_numpy()
        mcPosRsoErr[:,:,simNum] = mcLog[simNum].fswNavPosRsoEciErr.to_numpy()
        mcVelRsoErr[:,:,simNum] = mcLog[simNum].fswNavVelRsoEciErr.to_numpy()
        mcOeRsoErr[:,:,simNum] = mcLog[simNum].oeRso.to_numpy() - mcLog[simNum].fswNavOeRso.to_numpy()
        mcNrmInvnSqr[:,simNum] = mcLog[simNum].fswNavMeasNrmInvnSqr.to_numpy()
        mcNrmEstErrSqr[:,simNum] = mcLog[simNum].fswNavFltrNrmEstErrSqr.to_numpy()
    # Wrap angles to [-pi, pi)
    mcDoeErr[2:6,:,:]   = (mcDoeErr[2:6,:,:]   + np.pi) % (2 * np.pi) - np.pi
    mcOeVehErr[2:6,:,:] = (mcOeVehErr[2:6,:,:] + np.pi) % (2 * np.pi) - np.pi
    mcOeRsoErr[2:6,:,:] = (mcOeRsoErr[2:6,:,:] + np.pi) % (2 * np.pi) - np.pi

    print("MC: Computing Statistics ...")
    # Error standard deviations
    relPosRectRicErrStd = np.std(mcRelPosRectRicErr, axis = 2)
    relVelRectRicErrStd = np.std(mcRelVelRectRicErr, axis = 2)
    relPosCurvRicErrStd = np.std(mcRelPosCurvRicErr, axis = 2)
    relVelCurvRicErrStd = np.std(mcRelVelCurvRicErr, axis = 2)
    doeErrStd = np.std(mcDoeErr, axis = 2)
    posVehErrStd = np.std(mcPosVehErr, axis = 2)
    velVehErrStd = np.std(mcVelVehErr, axis = 2)
    oeVehErrStd = np.std(mcOeVehErr, axis = 2)
    posRsoErrStd = np.std(mcPosRsoErr, axis = 2)
    velRsoErrStd = np.std(mcVelRsoErr, axis = 2)
    oeRsoErrStd = np.std(mcOeRsoErr, axis = 2)
    # RMSE
    relPosRectRicRMSErr = rootMeanSum(mcRelPosRectRicErr)
    relVelRectRicRMSErr = rootMeanSum(mcRelVelRectRicErr)
    relPosCurvRicRMSErr = rootMeanSum(mcRelPosCurvRicErr)
    relVelCurvRicRMSErr = rootMeanSum(mcRelVelCurvRicErr)
    posVehRMSErr = rootMeanSum(mcPosVehErr)
    velVehRMSErr = rootMeanSum(mcVelVehErr)
    posRsoRMSErr = rootMeanSum(mcPosRsoErr)
    velRsoRMSErr = rootMeanSum(mcVelRsoErr)
    # Consistency
    avgNrmInvnSqr = np.mean(mcNrmInvnSqr, axis = 1)
    avgNrmEstErrSqr = np.mean(mcNrmEstErrSqr, axis = 1)

    print("MC: Packing Statistics ...")
    mcStats = {
        "time": mcLog[0].time.loc['simTime',:].to_numpy(),
        "relPosRectRicErrStd": relPosRectRicErrStd,
        "relVelRectRicErrStd": relVelRectRicErrStd,
        "relPosCurvRicErrStd": relPosCurvRicErrStd,
        "relVelCurvRicErrStd": relVelCurvRicErrStd,
        "doeErrStd": doeErrStd,
        "posVehErrStd": posVehErrStd,
        "velVehErrStd": velVehErrStd,
        "oeVehErrStd": oeVehErrStd,
        "posRsoErrStd": posRsoErrStd,
        "velRsoErrStd": velRsoErrStd,
        "oeRsoErrStd": oeRsoErrStd,
        "relPosRectRicRMSErr": relPosRectRicRMSErr,
        "relVelRectRicRMSErr": relVelRectRicRMSErr,
        "relPosCurvRicRMSErr": relPosCurvRicRMSErr,
        "relVelCurvRicRMSErr": relVelCurvRicRMSErr,
        "posVehRMSErr": posVehRMSErr,
        "velVehRMSErr": velVehRMSErr,
        "posRsoRMSErr": posRsoRMSErr,
        "velRsoRMSErr": velRsoRMSErr,
        "avgNrmInvnSqr": avgNrmInvnSqr,
        "avgNrmEstErrSqr": avgNrmEstErrSqr
    }
    
    return mcStats

def rootMeanSum(mcErr):
    return np.sqrt(np.mean(np.sum(mcErr**2, axis = 0), axis = 1))