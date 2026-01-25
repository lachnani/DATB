# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 22:21:40 2025

@author: Hakim Lachnani
"""

import os

import yaml
import tkinter as tk
from tkinter import filedialog
import ntpath
import pickle
import numpy as np

from dynamics import orbit as orb

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

file_path = filedialog.askopenfilename(initialdir = os.getcwd())
if ntpath.basename(file_path)[-4:] == '.pkl':
    """ We are loading a log file, not running a fresh sim """
    log_file_path = file_path
    with open(log_file_path, 'rb') as inp:
        log = pickle.load(inp)
        initConds = log[0]
        opsaTruth = log[1]
        tKovBreachTruth = log[2]
        opsaProp = log[3]
        opsaPropSpeed = log[4]
        opsaPropSpeedAvg = log[5]
        opsaPropErr = log[6]
        opsaPropErrAvg = log[7]
        opsaKorim = log[8]
        opsaKorimSpeed = log[9]
        opsaKorimSpeedAvg = log[10]
        opsaKorimErr = log[11]
        opsaKorimErrAvg = log[12]
        summaryTable = log[13]
    