# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:12:19 2026

@author: Hakim Lachnani
"""

import numpy as np

def computeMcStats(mcLog, numSims):
    print("MC: Accumulating Data ...")
    sim_len = len(mcLog[0].simTime)
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
    mcRelCovDet = np.zeros((2,sim_len,numSims))
    mcRelCovTrace = np.zeros((2,sim_len,numSims))
    mcRelCovEig = np.zeros((6,sim_len,numSims))
    mcRelCovFro = np.zeros((2,sim_len,numSims))
    mcVehCovDet = np.zeros((2,sim_len,numSims))
    mcVehCovTrace = np.zeros((2,sim_len,numSims))
    mcVehCovEig = np.zeros((6,sim_len,numSims))
    mcVehCovFro = np.zeros((2,sim_len,numSims))
    mcRsoCovDet = np.zeros((2,sim_len,numSims))
    mcRsoCovTrace = np.zeros((2,sim_len,numSims))
    mcRsoCovEig = np.zeros((6,sim_len,numSims))
    mcRsoCovFro = np.zeros((2,sim_len,numSims))
    mcRelNrmEstErrSqr = np.zeros((sim_len,numSims))
    mcVehNrmEstErrSqr = np.zeros((sim_len,numSims))
    mcRsoNrmEstErrSqr = np.zeros((sim_len,numSims))
    for simNum in range(numSims):
        mcRelPosRectRicErr[:,:,simNum] = mcLog[simNum].fswNavRelPosRectRicErr
        mcRelVelRectRicErr[:,:,simNum] = mcLog[simNum].fswNavRelVelRectRicErr
        mcRelPosCurvRicErr[:,:,simNum] = mcLog[simNum].fswNavRelPosCurvRicErr
        mcRelVelCurvRicErr[:,:,simNum] = mcLog[simNum].fswNavRelVelCurvRicErr
        mcDoeErr[:,:,simNum] = mcLog[simNum].fswNavDoeErr
        mcPosVehErr[:,:,simNum] = mcLog[simNum].fswNavPosVehEciErr
        mcVelVehErr[:,:,simNum] = mcLog[simNum].fswNavVelVehEciErr
        mcOeVehErr[:,:,simNum] = mcLog[simNum].oeVeh - mcLog[simNum].fswNavOeVeh
        mcPosRsoErr[:,:,simNum] = mcLog[simNum].fswNavPosRsoEciErr
        mcVelRsoErr[:,:,simNum] = mcLog[simNum].fswNavVelRsoEciErr
        mcOeRsoErr[:,:,simNum] = mcLog[simNum].oeRso - mcLog[simNum].fswNavOeRso
        mcRelCovDet[:,:,simNum] = mcLog[simNum].fswNavCovRelRicDet
        mcRelCovTrace[:,:,simNum] = mcLog[simNum].fswNavCovRelRicTrace
        mcRelCovEig[:,:,simNum] = mcLog[simNum].fswNavCovRelRicEig
        mcRelCovFro[:,:,simNum] = mcLog[simNum].fswNavCovRelRicFro
        mcVehCovDet[:,:,simNum] = mcLog[simNum].fswNavCovVehEciDet
        mcVehCovTrace[:,:,simNum] = mcLog[simNum].fswNavCovVehEciTrace
        mcVehCovEig[:,:,simNum] = mcLog[simNum].fswNavCovVehEciEig
        mcVehCovFro[:,:,simNum] = mcLog[simNum].fswNavCovVehEciFro
        mcRsoCovDet[:,:,simNum] = mcLog[simNum].fswNavCovRsoEciDet
        mcRsoCovTrace[:,:,simNum] = mcLog[simNum].fswNavCovRsoEciTrace
        mcRsoCovEig[:,:,simNum] = mcLog[simNum].fswNavCovRsoEciEig
        mcRsoCovFro[:,:,simNum] = mcLog[simNum].fswNavCovRsoEciFro
        mcNrmInvnSqr[:,simNum] = mcLog[simNum].fswNavMeasNrmInvnSqr
        mcRelNrmEstErrSqr[:,simNum] = mcLog[simNum].fswNavRelNrmEstErrSqr
        mcVehNrmEstErrSqr[:,simNum] = mcLog[simNum].fswNavVehNrmEstErrSqr
        mcRsoNrmEstErrSqr[:,simNum] = mcLog[simNum].fswNavRsoNrmEstErrSqr
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
    # Covariance
    avgRelCovDet = np.mean(mcRelCovDet, axis = 2)
    avgRelCovTrace = np.mean(mcRelCovTrace, axis = 2)
    avgRelCovEig = np.mean(mcRelCovEig, axis = 2)
    avgRelCovFro = np.mean(mcRelCovFro, axis = 2)
    avgVehCovDet = np.mean(mcVehCovDet, axis = 2)
    avgVehCovTrace = np.mean(mcVehCovTrace, axis = 2)
    avgVehCovEig = np.mean(mcVehCovEig, axis = 2)
    avgVehCovFro = np.mean(mcVehCovFro, axis = 2)
    avgRsoCovDet = np.mean(mcRsoCovDet, axis = 2)
    avgRsoCovTrace = np.mean(mcRsoCovTrace, axis = 2)
    avgRsoCovEig = np.mean(mcRsoCovEig, axis = 2)
    avgRsoCovFro = np.mean(mcRsoCovFro, axis = 2)
    # Consistency
    avgNrmInvnSqr = np.mean(mcNrmInvnSqr, axis = 1)
    avgRelNrmEstErrSqr = np.mean(mcRelNrmEstErrSqr, axis = 1)
    avgVehNrmEstErrSqr = np.mean(mcVehNrmEstErrSqr, axis = 1)
    avgRsoNrmEstErrSqr = np.mean(mcRsoNrmEstErrSqr, axis = 1)

    print("MC: Packing Statistics ...")
    mcStats = {
        "time": mcLog[0].simTime,
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
        "avgRelCovDet": avgRelCovDet,
        "avgRelCovTrace": avgRelCovTrace,
        "avgRelCovEig": avgRelCovEig,
        "avgRelCovFro": avgRelCovFro,
        "avgVehCovDet": avgVehCovDet,
        "avgVehCovTrace": avgVehCovTrace,
        "avgVehCovEig": avgVehCovEig,
        "avgVehCovFro": avgVehCovFro,
        "avgRsoCovDet": avgRsoCovDet,
        "avgRsoCovTrace": avgRsoCovTrace,
        "avgRsoCovEig": avgRsoCovEig,
        "avgRsoCovFro": avgRsoCovFro,
        "avgNrmInvnSqr": avgNrmInvnSqr,
        "avgRelNrmEstErrSqr": avgRelNrmEstErrSqr,
        "avgVehNrmEstErrSqr": avgVehNrmEstErrSqr,
        "avgRsoNrmEstErrSqr": avgRsoNrmEstErrSqr
    }
    
    return mcStats

def rootMeanSum(mcErr):
    return np.sqrt(np.mean(np.sum(mcErr**2, axis = 0), axis = 1))