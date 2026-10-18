# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 18:53:07 2025

@author: Hakim Lachnani
"""

import os

import tkinter as tk
from tkinter import filedialog
import importlib
import ntpath
import time
from datetime import datetime
import pickle
import numpy as np
import inspect
import copy
from pathlib import Path

import simulator
import parser
from plotting import plotMonteCarlo as pltMC

        
""" 
Runs Kalman Filter Monte Carlo Simulations

"""

print("-----------------------------------------------------")
print("-----------------------------------------------------")
print("D.A.T.B: Differential Astrodynamics Test Bench")
print("October 2026")
print("-----------------------------------------------------")

root = tk.Tk()
root.withdraw()

""" Load the config parameters """
print("Choose Monte Carlo yaml file")
cfg, cfg_file = parser.loadFile('Monte Carlo Configuration')
        
numSims = cfg["numSims"]
settings = cfg["sim"]
if settings["fsw"]["status"] == False:
    raise ValueError("Settings file must have FSW enabled")
settings["fsw"]["initErr"] = True # Initialize error
formation = cfg["init"]
frm = parser.parseFormation(formation)
fsw_config = cfg["fsw"]
fsw = parser.parseFlightSoftwareInit(formation['epoch'], settings['fsw'])

""" Load the scenario """
print("Choose scenario file")
scr_file_path = filedialog.askopenfilename(defaultextension = '.py',
                                       initialdir = os.getcwd())
scr_file = ntpath.basename(scr_file_path)[:-3]
scr = importlib.import_module('scenarios.' + scr_file)
scrArguments = inspect.getfullargspec(scr.scenario).args

if 'wptTbl' in scrArguments:
    """ Generate the Waypoint Table """   
    wptGen, wptGen_file = parser.loadFile('waypoint table generator')
    wptTbl, frm = parser.parseWptTbl(wptGen, frm)
    if settings["simDuration"] < wptTbl.t[-1]:
        settings["simDuration"] = wptTbl.t[-1]

# Define logs
mcLog = []

# Run the sims
print("MC: Running Sims ...")
start_time = time.time()
for simNum in range(numSims):
    print("MC: Running Sim {}".format(simNum+1))
    # Copy and initialize structures
    frm_case = copy.deepcopy(frm)
    fsw_case = copy.deepcopy(fsw)
    fsw_case = parser.parseFlightSoftware(frm_case, fsw_case, fsw_config, True)
    
    # Initialize the sim
    sim_case = simulator.Simulator(settings, frm_case, fsw_case, quiet = True)
    
    # Run the Sim
    if 'wptTbl' in scrArguments:
        scr.scenario(sim_case, wptTbl)
    else:
        scr.scenario(sim_case)
    
    # Save relevant log information
    mcLog.append(sim_case.log)
    
    
end_time = time.time()
elapsed_time = end_time - start_time
print("MC: Simulations completed in %.1f seconds" % elapsed_time)
print("MC: Average Realtime factor: %.2f" % (settings["simDuration"]*numSims/elapsed_time))
print("-----------------------------------------------------")

""" Compute statistics """
print("MC: Accumulating Data ...")
sim_len = len(mcLog[0].time.loc['simTime',:])
mcRelPosRectRicErr = np.zeros((3,sim_len,numSims))
mcRelVelRectRicErr = np.zeros((3,sim_len,numSims))
mcRelPosCurvRicErr = np.zeros((3,sim_len,numSims))
mcRelVelCurvRicErr = np.zeros((3,sim_len,numSims))
mcOeVehErr = np.zeros((6,sim_len,numSims))
mcOeRsoErr = np.zeros((6,sim_len,numSims))
for simNum in range(numSims):
    mcRelPosRectRicErr[:,:,simNum] = mcLog[simNum].fswNavRelPosRectRicErr.to_numpy()
    mcRelVelRectRicErr[:,:,simNum] = mcLog[simNum].fswNavRelVelRectRicErr.to_numpy()
    mcRelPosCurvRicErr[:,:,simNum] = mcLog[simNum].fswNavRelPosCurvRicErr.to_numpy()
    mcRelVelCurvRicErr[:,:,simNum] = mcLog[simNum].fswNavRelVelCurvRicErr.to_numpy()
    mcOeVehErr[:,:,simNum] = mcLog[simNum].oeVeh.to_numpy() - mcLog[simNum].fswNavOeVeh.to_numpy()
    mcOeRsoErr[:,:,simNum] = mcLog[simNum].oeRso.to_numpy() - mcLog[simNum].fswNavOeRso.to_numpy()

print("MC: Computing Statistics ...")
relPosRectRicErrStd = np.std(mcRelPosRectRicErr, axis = 2)
relVelRectRicErrStd = np.std(mcRelVelRectRicErr, axis = 2)
relPosCurvRicErrStd = np.std(mcRelPosCurvRicErr, axis = 2)
relVelCurvRicErrStd = np.std(mcRelVelCurvRicErr, axis = 2)
oeVehErrStd = np.std(mcOeVehErr, axis = 2)
oeRsoErrStd = np.std(mcOeRsoErr, axis = 2)

print("MC: Packing Statistics ...")
mcStats = {
    "time": mcLog[0].time.loc['simTime',:].to_numpy(),
    "relPosRectRicErrStd": relPosRectRicErrStd,
    "relVelRectRicErrStd": relVelRectRicErrStd,
    "relPosCurvRicErrStd": relPosCurvRicErrStd,
    "relVelCurvRicErrStd": relVelCurvRicErrStd,
    "oeVehErrStd": oeVehErrStd,
    "oeRsoErrStd": oeRsoErrStd
}

""" Save Log """
print("MC: Saving Data ...")
log_time = datetime.now()
log_time_str = dts = log_time.strftime("%Y%h%d_%H%M")
log_folder = cfg_file + '_' \
           + scr_file + '_' \
           + log_time_str
log_base_path = os.path.join("log", log_folder)
Path(log_base_path).mkdir(parents=True, exist_ok=True)

log_file_path = os.path.join(log_base_path, 'log.pkl')
with open(log_file_path, 'wb') as outp:
    pickle.dump(mcLog, outp)

""" Plots """
print("MC: Plotting Data ...")
pltMC.plotAll(mcLog, mcStats, log_base_path, log_time_str, settings)
    


        

    