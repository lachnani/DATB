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
from copy import deepcopy

import simulator
import parser
from dynamics import formation
from planning import wptTbl as wt

        
""" 
Analyzes waypoint table plan. 

Simulates a waypoint plan in different configurations:
    1. Nominal plan simulation
    2. Missed burn branch simulations
    3. Burn dispersion Monte Carlo simulations
    4. Full FSW emulation Monte Carlo simulations

"""

print("-----------------------------------------------------")
print("-----------------------------------------------------")
print("D.A.T.B: Differential Astrodynamics Test Bench")
print("-----------------------------------------------------")
print("Waypoint Plan Analysis")
print("-----------------------------------------------------")

root = tk.Tk()
root.withdraw()

"""
Generate the Formation from the Initial State
"""

print("Choose initial state yaml file")
state_file_path = filedialog.askopenfilename(defaultextension = '.yaml',
                                       initialdir = os.getcwd())
state_file = ntpath.basename(state_file_path)[:-5]
with open(state_file_path) as stream:
    try:
        formationCfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

frm_base = parser.parseFormation(formationCfg)


"""
Generate the Waypoint Table
"""

print("Choose waypoint table generator yaml file")        
cfg_file_path = filedialog.askopenfilename(defaultextension = '.yaml',
                                       initialdir = os.getcwd())
cfg_file = ntpath.basename(cfg_file_path)[:-5]
with open(cfg_file_path) as stream:
    try:
        cfg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

generator = cfg.pop("generator")
if generator["type"] == "CLROE":
    wptTbl = wt.waypointTable()
    wptTbl.genClroeWpts(
        parser.unpackClroe(generator["L"]),
        generator["meanMotion"],
        cfg["tblStartTime"],
        cfg["tblEndTime"],
        cfg["numWpts"])
    
    if cfg["relStateNatural"]:
        # Replace the formation deputy position with the natural path.
        # Assume curvilinear for long-range accuracy
        frm_base = formation.Formation(
            frm_base.chief, 
            frm_base.deputy, 
            parser.unpackClroe(generator["L"]),
            "FORMATION_CHIEF_ANCHOR",
            "RELSTATE_CURV_CLROE",
            frm_base.chief.pert)
        
"""
Configure Simulation Settings
"""
        
print("Choose analysis settings yaml file")        
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
    "simDuration": wptTbl.t[-1],
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
frm_case = deepcopy(frm_base)

print("Choose scenario file")
scr_file_path = filedialog.askopenfilename(defaultextension = '.py',
                                       initialdir = os.getcwd())
scr_file = ntpath.basename(scr_file_path)[:-3]
scr = importlib.import_module('scenarios.' + scr_file)

# Initialize the simulation
sim = simulator.Simulator(settings_case, 
                          frm_case, 
                          0, 
                          quiet = True)

# Run the sim
scr.scenario(sim, wptTbl)


