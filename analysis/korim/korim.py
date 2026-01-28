# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 18:53:07 2025

@author: Hakim Lachnani
"""

import os

import yaml
import tkinter as tk
from tkinter import filedialog
import ntpath
import time
from datetime import datetime
import pickle
from pathlib import Path
import pandas as pd
import numpy as np
from numpy import random as rand
import matplotlib.pyplot as plt
from itertools import product, combinations

import simulator
import parser
from dynamics import orbit as orb
from planning import passiveSafety as ps

        
""" 
Runs KORIM Monte Carlo Scenarios

"""

print("-----------------------------------------------------")
print("-----------------------------------------------------")
print("D.A.T.B: Differential Astrodynamics Test Bench")
print("Version 0.4 --- July 2025")
print("-----------------------------------------------------")

root = tk.Tk()
root.withdraw()

print("Choose KORIM Monte Carlo yaml file")
""" Load the config parameters """
cfg_file_path = filedialog.askopenfilename(defaultextension = '.yaml',
                                       initialdir = os.getcwd())
cfg_file = ntpath.basename(cfg_file_path)[:-5]
with open(cfg_file_path) as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        
""" Create the default settings dictionary """
print("MC: Loading Dynamics ...")
dynamics = cfg["dynamics"]
log = {"status": True, "dt": 1.0}
visualizer = {"plots": False, "animations": False}
fsw = {"status": False, "dt": 1.0}
orbit = {
    "environments": True,
    "elements": True
    }
formation = {
    "relStates": True,
    "environments": False,
    "measurements": False,
    "orbit": orbit
    }
settings_base = {
    "name": "KORIM",
    "dynamics": cfg["dynamics"],
    "simDuration": 1,
    "log": log,
    "visualizer": visualizer,
    "fsw": fsw,
    "formation": formation
    }
cfg.pop('dynamics')
        
"""
Run the truth simulations

Define the remaining sim parameters according to the max oth and the chief 
semimajor axis:
    - log dt
    - simDuration
For each sim:
    - Compute the initial relative state
"""

kovLim = np.array([dynamics["kov"]["R"],dynamics["kov"]["I"],dynamics["kov"]["C"]])

# Default formation dictionary
print("MC: Loading base formation ...")
frm_base = {
    "epoch": cfg["epoch"],
    "chief": cfg["chief"],
    "deputy": 0,
    "relState":{
        "state":{
            "A0": 0.0,
            "alpha": 0.0,
            "xOff": 0.0,
            "yOff": 0.0,
            "B0": 0.0,
            "beta": 0.0
            },
        },
    "frmType": "FORMATION_CHIEF_ANCHOR",
    "relStateType": cfg["relStateType"],
    "pert": cfg["pert"]
    }

# Compute sim duration
oths = np.array(cfg["oth"].split(), dtype = 'f')
numOth = oths.size
orbPeriod = np.ceil(2*np.pi*np.sqrt(frm_base["chief"]["state"]["a"]**3/orb.MU_EARTH))
meanMotion = 2*np.pi/orbPeriod
simDuration = np.floor(np.max(oths))*orbPeriod

# Set sim settings
settings_case = settings_base.copy()
settings_case["simDuration"] = simDuration
settings_case["log"]["dt"] = orbPeriod

# Define logs
initConds = np.zeros((cfg["numSims"],6))
opsaTruth = np.ones((cfg["numSims"],numOth))
tKovBreachTruth = np.zeros((cfg["numSims"],))

# Run the sims
print("MC: Running Sims ...")
start_time = time.time()
for simNum in range(cfg["numSims"]):
    print("MC: Running Sim {}".format(simNum+1))
    # Compute initial CLROEs using uniform distribution
    frm_sim = frm_base.copy()
    frm_sim["relState"]["state"]["A0"] = rand.uniform(
        cfg["relState"]["limits"]["A0"]["min"],
        cfg["relState"]["limits"]["A0"]["max"])
    frm_sim["relState"]["state"]["alpha"] = rand.uniform(
        np.deg2rad(cfg["relState"]["limits"]["alpha"]["min"]),
        np.deg2rad(cfg["relState"]["limits"]["alpha"]["max"]))
    frm_sim["relState"]["state"]["xOff"] = rand.uniform(
        cfg["relState"]["limits"]["xOff"]["min"],
        cfg["relState"]["limits"]["xOff"]["max"])
    frm_sim["relState"]["state"]["yOff"] = rand.uniform(
        cfg["relState"]["limits"]["yOff"]["min"],
        cfg["relState"]["limits"]["yOff"]["max"])
    frm_sim["relState"]["state"]["B0"] = rand.uniform(
        cfg["relState"]["limits"]["B0"]["min"],
        cfg["relState"]["limits"]["B0"]["max"])
    frm_sim["relState"]["state"]["beta"] = rand.uniform(
        np.deg2rad(cfg["relState"]["limits"]["beta"]["min"]),
        np.deg2rad(cfg["relState"]["limits"]["beta"]["max"]))
    
    # Initialize the sim
    sim = simulator.Simulator(settings_case, parser.parseFormation(frm_sim), 0, quiet = True)
    if (frm_base["relStateType"] == "RELSTATE_CURV_CLROE"):
        initConds[simNum,:] = sim.frm.curvClroe
    elif (frm_base["relStateType"] == "RELSTATE_RECT_CLROE"):
        initConds[simNum,:] = sim.frm.rectClroe
    
    # Run the Sim
    sim.run(simDuration)
    
    # Log OPSA truth data
    tKovBreachTruth[simNum] = sim.log.tKovBreach
    if not sim.log.kovBreach:
        opsaTruth[simNum,:] = np.ones((1,numOth))
    else:
        for othNum in range(numOth):     
            oth = oths[othNum]
            if sim.log.tKovBreach < oth*orbPeriod:
                opsaTruth[simNum,othNum] = 0
    
    
end_time = time.time()
elapsed_time = end_time - start_time
print("MC: Simulations completed in %.1f seconds" % elapsed_time)
print("MC: Average Realtime factor: %.2f" % (simDuration*cfg["numSims"]/elapsed_time))
print("-----------------------------------------------------")
    
# Run and Clock Propagation
print("MC: Running Propagation OPSA ...")
propInts = np.array(cfg["propInt"].split(), dtype = 'f')
numPropInt = propInts.size
opsaProp = np.zeros((numPropInt,cfg["numSims"],numOth))
opsaPropSpeedAvg = np.zeros((numPropInt,numOth))
opsaPropErr = np.zeros((numPropInt,cfg["numSims"],numOth))
opsaPropErrAvg = np.zeros((numPropInt,numOth))
for propIntNum in range(numPropInt):
    propInt = propInts[propIntNum]
    for othNum in range(numOth):
        oth = oths[othNum]  
        start_time_avg = time.time()
        for simNum in range(cfg["numSims"]):
            opsaProp[propIntNum, simNum, othNum] = ps.clroeProp(initConds[simNum,:],
                                                                          meanMotion,
                                                                          oth*orbPeriod,
                                                                          propInt,
                                                                          kovLim)
        end_time_avg = time.time()
        elapsed_time = end_time_avg - start_time_avg
        opsaPropSpeedAvg[propIntNum, othNum] = elapsed_time/cfg["numSims"]
        print(f"MC: Average Prop OPSA Speed over {oth:.1f} period OTH with {propInt/60:.1f} minute interval: {opsaPropSpeedAvg[propIntNum, othNum]:.2E}")
        opsaPropErr[propIntNum, :, othNum] = opsaProp[propIntNum, :, othNum] - opsaTruth[:, othNum]
        opsaPropErrAvg[propIntNum,othNum] = np.count_nonzero(opsaPropErr[propIntNum, :, othNum]) / cfg["numSims"]
        print(f"MC: Average Prop OPSA Accuracy over {oth:.1f} period OTH with {propInt/60:.1f} minute interval: {100*(1-opsaPropErrAvg[propIntNum,othNum]):.2f}%")
print("-----------------------------------------------------")
        
# Run and Clock KORIM
print("MC: Running KORIM OPSA ...")
opsaKorim = np.zeros((cfg["numSims"],numOth))
opsaKorimSpeedAvg = np.zeros((numOth,))
opsaKorimErr = np.zeros((cfg["numSims"],numOth))
opsaKorimErrAvg = np.zeros((numOth,))
for othNum in range(numOth):
    oth = oths[othNum] 
    start_time_avg = time.time()
    for simNum in range(cfg["numSims"]):
        opsaKorim[simNum,othNum] = ps.korim(initConds[simNum,:],
                                                      oth,
                                                      kovLim)

    end_time_avg = time.time()
    elapsed_time = end_time_avg - start_time_avg
    opsaKorimSpeedAvg[othNum] = elapsed_time/cfg["numSims"]
    print(f"MC: Average KORIM OPSA Speed over {oth:.1f} period OTH: {opsaKorimSpeedAvg[othNum]:.2E}")
    opsaKorimErr[:,othNum] = opsaKorim[:,othNum] - opsaTruth[:,othNum]
    opsaKorimErrAvg[othNum] = np.count_nonzero(opsaKorimErr[:,othNum]) / cfg["numSims"]
    print(f"MC: Average Prop OPSA Accuracy over {oth:.1f} period OTH: {100*(1-opsaKorimErrAvg[othNum]):.2f}%")
print("-----------------------------------------------------")

# Construct final case table
print("MC: Generating Summary Table")
summaryTable = pd.DataFrame({'Korim Runtime': opsaKorimSpeedAvg,
                   'Korim Accuracy': 100*(1-opsaKorimErrAvg)},
                  index=cfg["oth"].split())
for propIntNum in range(numPropInt):
    propInt = propInts[propIntNum]
    propStr = f"{propInt/60:.0f}m Prop Runtime"
    summaryTable.insert(2*propIntNum + 2, propStr, opsaPropSpeedAvg[propIntNum, :])
    propStr = f"{propInt/60:.0f}m Prop Accuracy"
    summaryTable.insert(2*propIntNum + 3, propStr, 100*(1-opsaPropErrAvg[propIntNum, :]))
print(summaryTable)
print("-----------------------------------------------------")

# Save relevant data
print("MC: Saving Monte Carlo data log")
log_time = datetime.now()
log_time_str = dts = log_time.strftime("%Y%h%d_%H%M")
log_folder = 'korim_' + cfg_file
log_file = "log_" + log_time_str
log_inner_path = os.path.join("log", log_folder)
log_file_path = os.path.join(log_inner_path, log_file + '.' + 'pkl')
Path(log_inner_path).mkdir(parents=True, exist_ok=True)
with open(log_file_path, 'wb') as outp:
    pickle.dump([initConds,
                 opsaTruth,
                 tKovBreachTruth,
                 opsaProp,
                 opsaPropSpeedAvg,
                 opsaPropErr,
                 opsaPropErrAvg,
                 opsaKorim,
                 opsaKorimSpeedAvg,
                 opsaKorimErr,
                 opsaKorimErrAvg,
                 summaryTable], outp)
print("-----------------------------------------------------")

# For the cases where KORIM fails, we want to show that the true trajectory 
# barely clips the edges of the KOV
print("MC: Checking False-Positive results ...")
opsaKorimFP = np.argwhere(opsaKorimErr > 0) 
unique, unique_indices = np.unique(opsaKorimFP[:,0], return_index = True)  
opsaKorimFP = opsaKorimFP[unique_indices,:]
numFP = np.shape(opsaKorimFP)[0]
opsaKorimFPMinRange = np.zeros((numFP,))
settings_fp = settings_case.copy()
settings_fp['log']['dt'] = 30
for othNum in range(numOth):
    oth = oths[othNum] 
    # Initialize plot
    korimFig = plt.figure()
    # plt.title(f'KORIM False-Positive Trajectories: {oth:.0f} Orbit OTH')
    ax = korimFig.add_subplot(projection="3d")
    origin = ax.scatter(0, 0, 0, color='b', linewidths = 3)
    ax.set_xlabel("Curvilinear In-Track [km]")
    ax.set_ylabel("Curvilinear Radial [km]")
    ax.set_zlabel("Curvilinear Cross-Track [km]")
    # Create cube
    r = [-kovLim[0], kovLim[0]]
    i = [-kovLim[1], kovLim[1]]
    c = [-kovLim[2], kovLim[2]]
    for s, e in combinations(np.array(list(product(i, r, c))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0] or np.sum(np.abs(s-e)) == i[1]-i[0] or np.sum(np.abs(s-e)) == c[1]-c[0]:
            ax.plot3D(*zip(s, e), color="b")
    for i in range(numFP):
        if othNum == opsaKorimFP[i,1]:
            simNum = opsaKorimFP[i,0]
            settings_fp["simDuration"] = oths[opsaKorimFP[i,1]]*orbPeriod
            # Rerun each false positive simulation with a 30s log rate
            print("MC: Re-running Sim {}".format(simNum+1))
            # Compute initial CLROEs using uniform distribution
            frm_sim = frm_base.copy()
            frm_sim["relState"]["state"]["A0"] = initConds[simNum,0]
            frm_sim["relState"]["state"]["alpha"] = initConds[simNum,1]
            frm_sim["relState"]["state"]["xOff"] = initConds[simNum,2]
            frm_sim["relState"]["state"]["yOff"] = initConds[simNum,3]
            frm_sim["relState"]["state"]["B0"] = initConds[simNum,4]
            frm_sim["relState"]["state"]["beta"] = initConds[simNum,5]
            
            # Initialize the sim
            sim = simulator.Simulator(settings_fp, parser.parseFormation(frm_sim), 0, quiet = True)
            
            # Run the Sim
            sim.run(simDuration)
            
            # Plot the result
            iPlot = sim.log.relPosCurvRic.loc['I',0:sim.log.i].to_numpy()
            iBreachInd = np.abs(iPlot) <= kovLim[1]
            rPlot = sim.log.relPosCurvRic.loc['R',0:sim.log.i].to_numpy()
            rBreachInd = np.abs(rPlot) <= kovLim[0]
            cPlot = sim.log.relPosCurvRic.loc['C',0:sim.log.i].to_numpy()
            cBreachInd = np.abs(cPlot) <= kovLim[2]
            
            allBreachInd = np.bitwise_and(np.bitwise_and(iBreachInd, rBreachInd), cBreachInd)
            
            iPlot[~allBreachInd] = None
            rPlot[~allBreachInd] = None
            cPlot[~allBreachInd] = None
            
            if any(allBreachInd):
                ax.plot(iPlot[allBreachInd], 
                        rPlot[allBreachInd], 
                        cPlot[allBreachInd],
                        linewidth=1.0, color = 'r')[0]
                
                # Save the minimum range
                opsaKorimFPMinRange[i] = np.min(np.linalg.norm(np.c_[rPlot[allBreachInd],iPlot[allBreachInd],cPlot[allBreachInd]], axis=1))
            else:
                opsaKorimFPMinRange[i] = np.max(kovLim)

    # Clean up the plot
    plt.grid()
    plt.gca().invert_xaxis()
    ax.set_aspect('equal', adjustable = 'datalim')
        

    