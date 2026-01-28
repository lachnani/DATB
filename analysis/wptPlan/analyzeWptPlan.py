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
import importlib
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
Analyzes waypoint table plan. 

Simulates a waypoint plan in different configurations:
    1. Nominal plan simulation
    2. Missed burn branch simulations
    3. Burn dispersion Monte Carlo simulations
    4. Full FSW emulation Monte Carlo simulations

"""

def loadPlan():
    print("Choose waypoint plan .pkl file")
    file_path = filedialog.askopenfilename(initialdir = os.getcwd())
    if ntpath.basename(file_path)[-4:] == '.pkl':
        """ We are loading a log file, not running a fresh sim """
        log_file_path = file_path
        with open(log_file_path, 'rb') as inp:
            plan = pickle.load(inp)
            frm = plan[0]
            wptTbl = plan[1]
        return frm, wptTbl

print("-----------------------------------------------------")
print("-----------------------------------------------------")
print("D.A.T.B: Differential Astrodynamics Test Bench")
print("-----------------------------------------------------")
print("Waypoint Plan Analysis")
print("-----------------------------------------------------")

root = tk.Tk()
root.withdraw()

"""Load the waypoint plan"""
frm_base, wptTbl = loadPlan()
        
print("Choose analysis settings yaml file")        
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
    "environments": True,
    "measurements": True,
    "orbit": orbit
    }
settings_base = {
    "name": "analyzeWptPlan",
    "dynamics": cfg["dynamics"],
    "simDuration": 1,
    "log": log,
    "visualizer": visualizer,
    "fsw": fsw,
    "formation": formation
    }
cfg.pop('dynamics')

"""
Nominal plan simulation
"""

settings_case = settings_base.copy()
frm_case = frm_base.copy()

print("Choose scenario file")
scr_file_path = filedialog.askopenfilename(defaultextension = '.py',
                                       initialdir = os.getcwd())
scr_file = ntpath.basename(scr_file_path)[:-3]
scr = importlib.import_module('scenarios.' + scr_file)

# Initialize the simulation
sim = simulator.Simulator(settings_case, 
                          parser.parseFormation(frm_case), 
                          0, 
                          quiet = True)

# Run the sim
scr.scenario(sim, wptTbl)


