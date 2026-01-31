# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:38:52 2025

@author: Hakim Lachnani
"""

import os

import tkinter as tk
from tkinter import filedialog
import importlib
import ntpath
from datetime import datetime
import pickle
from pathlib import Path
import inspect

import simulator
import parser
import visualizer as vis

        
""" 
Running a Simulation:
    Load the configuration parameters:
        Simulator configuration parameters
        Initial state
        FSW paramaters
    Pick a scenario
    Run the sim

"""

print("-----------------------------------------------------")
print("-----------------------------------------------------")
print("D.A.T.B: Differential Astrodynamics Test Bench")
print("Version 0.4 --- January 2026")
print("-----------------------------------------------------")

root = tk.Tk()
root.withdraw()

""" Load the sim config parameters """
settings, sim_file = parser.loadFile('simulation parameter')
        
""" Load the initial state """
formation, state_file = parser.loadFile('initial state')

        
""" Load the FSW """
if settings["fsw"]["status"] ==  True:
    fsw_config, fsw_file = parser.loadFile('FSW')
else:
    fsw_config = 0
    fsw_file = 'noFSW'
    
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
    wptTbl, formation = parser.parseWptTbl(wptGen, formation)
    

""" Initialize the sim """
print("-----------------------------------------------------")
sim = simulator.Simulator(
    settings, 
    parser.parseFormation(formation), 
    fsw_config)


""" Run the sim """
scr.scenario(sim)


""" Log the data """
if settings["log"]["status"] ==  True:
    print("-----------------------------------------------------")
    print("Saving simulation data log")
    log_time = datetime.now()
    log_time_str = dts = log_time.strftime("%Y%h%d_%H%M")
    log_folder = sim_file + '_' \
               + state_file + '_' \
               + fsw_file + '_' \
               + scr_file 
    log_file = "log_" + log_time_str
    log_inner_path = os.path.join("log", log_folder)
    log_file_path = os.path.join(log_inner_path, log_file + '.' + 'pkl')
    Path(log_inner_path).mkdir(parents=True, exist_ok=True)
    with open(log_file_path, 'wb') as outp:
        pickle.dump(sim.log, outp)
        
    """ Visualize the data """
    if settings["visualizer"]["plots"] ==  True:
        print("Saving simulation plots")
        vis.plotAll(sim.log, log_inner_path, log_time_str, sim.settings)
        
    if settings["visualizer"]["animations"] ==  True:
        print("Saving simulation animations")
        vis.animAll(sim.log, log_inner_path, log_time_str, sim.settings)