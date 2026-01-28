# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 13:03:36 2026

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

from planning import wptTbl as wt
import parser

def savePlan(frm, wptTbl, path = None, name = None):
    print("Saving waypoint plan .pkl file")
    if name == None:
        nameTime = datetime.now()
        name = nameTime.strftime("%Y%h%d_%H%M")
    file = "plan_" + name
    if path == None:
        file_path = os.path.join(file + '.' + 'pkl')
    else: 
        file_path = os.path.join(path + "\\" + file + '.' + 'pkl')
    with open(file_path, 'wb') as plan:
        pickle.dump([frm, wptTbl], plan)
        
        
print("-----------------------------------------------------")
print("-----------------------------------------------------")
print("D.A.T.B: Differential Astrodynamics Test Bench")
print("-----------------------------------------------------")
print("Waypoint Plan Generation")
print("-----------------------------------------------------")

"""
Generate the Formation
"""


"""
Generate the Waypoint Table
"""

print("Choose waypoint table generator yaml file")        
""" Load the config parameters """
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
    wptTbl = wt.wptTbl()
    wptTbl.genClroeWpts(
        parser.unpackClroe(generator["L"]),
        generator["meanMotion"],
        cfg["tblStartTime"],
        cfg["tblEndTime"],
        cfg["numWPts"])
    
    
