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
from navigation import estimator
import kfMcUtils
from planning import dvTbl

""" 
Runs Kalman Filter Monte Carlo Comparison Simulations

"""

"""
Hard-coded inputs
"""

testFilterFormulations = [
    "DualInertial",
    "InertialRelative",
    "ChiefDecoupledDualInertial",
    "RelativeDecoupledInertialRelative"
    ]

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
        
if 'dvTbl' in scrArguments:
    dvt = dvTbl.deltaVTable()
    dvt.appendDvs(np.array(cfg['dvTbl']['t']),np.array(cfg['dvTbl']['dv']).reshape(3, 1),cfg['dvTbl']['frame'])
    acc = cfg['acc']
        
# Generate random input
print("MC: Generating Inputs ...")
# Initial Covariance
deputyCovInit = estimator.initCovFromRic(
    np.array([fsw_config['navigation']['deputyCov']['R'], fsw_config['navigation']['deputyCov']['I'], fsw_config['navigation']['deputyCov']['C'],
              fsw_config['navigation']['deputyCov']['Rv'],fsw_config['navigation']['deputyCov']['Iv'],fsw_config['navigation']['deputyCov']['Cv']]), 
    frm.deputy.r, 
    frm.deputy.v)
chiefCovInit = estimator.initCovFromRic(
    np.array([fsw_config['navigation']['chiefCov']['R'], fsw_config['navigation']['chiefCov']['I'], fsw_config['navigation']['chiefCov']['C'],
              fsw_config['navigation']['chiefCov']['Rv'],fsw_config['navigation']['chiefCov']['Iv'],fsw_config['navigation']['chiefCov']['Cv']]), 
    frm.chief.r, 
    frm.chief.v)
mcDeputyStateErr = np.random.multivariate_normal(np.zeros(6,),deputyCovInit, numSims).T
mcChiefStateErr = np.random.multivariate_normal(np.zeros(6,),chiefCovInit, numSims).T

# Run the sims
print("MC: Running Sims ...")
start_time = time.time()
mcStats = []
for fltrFrm in testFilterFormulations:
    print(f"MC: Running Formulation {fltrFrm}")
    # Overide filter formulation
    fsw_config_case = copy.deepcopy(fsw_config)
    fsw_config_case['navigation']['filterType'] = fltrFrm
    
    # Clear log
    mcLog = []
    
    # Run sims
    for simNum in range(numSims):
        print("MC: Running Sim {}".format(simNum+1))
        # Copy and initialize structures
        frm_case = copy.deepcopy(frm)
        fsw_case = copy.deepcopy(fsw)
        fsw_case = parser.parseFlightSoftware(frm_case, fsw_case, fsw_config_case, 
                                              initErr = False,
                                              deputyStateErr = mcDeputyStateErr[:,simNum],
                                              chiefStateErr = mcChiefStateErr[:,simNum])
        
        # Initialize the sim
        sim_case = simulator.Simulator(settings, frm_case, fsw_case, quiet = True)
        
        # Run the Sim
        if 'wptTbl' in scrArguments:
            scr.scenario(sim_case, wptTbl)
        elif 'dvTbl' in scrArguments:
            scr.scenario(sim_case, dvt, acc)
        else:
            scr.scenario(sim_case)
        
        # Save relevant log information
        mcLog.append(sim_case.log)
        
    # Compute statistics
    mcStats.append(kfMcUtils.computeMcStats(mcLog, numSims))
    
    
end_time = time.time()
elapsed_time = end_time - start_time
print("MC: Simulations completed in %.1f seconds" % elapsed_time)
print("MC: Average Realtime factor: %.2f" % (settings["simDuration"]*numSims*len(testFilterFormulations)/elapsed_time))
print("-----------------------------------------------------")


""" Save Log """
print("MC: Saving Data ...")
log_time = datetime.now()
log_time_str = dts = log_time.strftime("%Y%h%d_%H%M")
log_folder = 'kfComp' + '_' \
           + cfg_file + '_' \
           + scr_file + '_' \
           + log_time_str
log_base_path = os.path.join("log", log_folder)
Path(log_base_path).mkdir(parents=True, exist_ok=True)

log_file_path = os.path.join(log_base_path, 'log.pkl')
with open(log_file_path, 'wb') as outp:
    pickle.dump({'cfg': cfg,
                 'mcStats': mcStats},
                outp)

""" Plots """
print("MC: Plotting Data ...")
pltMC.plotComp(mcStats, testFilterFormulations, log_base_path, log_time_str)
    


        

    