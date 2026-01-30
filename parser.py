# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 18:33:02 2025

@author: Hakim Lachnani
"""

import os

import yaml
from tkinter import filedialog
import ntpath
import numpy as np

from dynamics import formation as frm
from dynamics import orbit as orb
from planning import wptTbl as wt

def loadFile(promptStr = ''):
    print("Choose " + promptStr + " yaml file")
    file_path = filedialog.askopenfilename(defaultextension = '.yaml',
                                           initialdir = os.getcwd())
    file = ntpath.basename(file_path)[:-5]
    with open(file_path) as stream:
        try:
            return yaml.safe_load(stream), file
        except yaml.YAMLError as exc:
            print(exc)

def unpackClroe(stateDict):
    return np.array([
        stateDict["A0"],
        np.deg2rad(stateDict["alpha"]),
        stateDict["xOff"],
        stateDict["yOff"],
        stateDict["B0"],
        np.deg2rad(stateDict["beta"])])

def unpackKepEl(stateDict):
    return np.array([
        stateDict["a"],
        stateDict["e"],
        np.deg2rad(stateDict["i"]),
        np.deg2rad(stateDict["RAAN"]),
        np.deg2rad(stateDict["argP"]),
        np.deg2rad(stateDict["M"])])

def unpackDoe(stateDict):
    return np.array([
        stateDict["da"],
        stateDict["de"],
        np.deg2rad(stateDict["di"]),
        np.deg2rad(stateDict["dRAAN"]),
        np.deg2rad(stateDict["dargP"]),
        np.deg2rad(stateDict["dM"])])

def unpackRelPosVel(stateDict):
    return np.array([
        stateDict["relPosR"],
        stateDict["relPosI"],
        stateDict["relPosC"],
        stateDict["relVelR"],
        stateDict["relVelI"],
        stateDict["relVelC"]])

def parseFormation(yaml):
    """
    Parse Formation yaml file

    Parameters
    ----------
    yaml : dictionary from yaml

    Returns
    -------
    formation object

    """
    
    chief = 0
    deputy = 0
    relState = 0

    # Define the anchor
    if ((yaml["frmType"] == "FORMATION_CHIEF_ANCHOR") or
        (yaml["frmType"] == "FORMATION_DUAL_ABS")):
        # Define the chief
        if yaml["chief"]["stateType"] == "STATE_KEPEL":
            chief = orb.Orbit(
                yaml["epoch"], 
                unpackKepEl(yaml["chief"]["state"]),
                stateType = yaml["chief"]["stateType"],
                pert = yaml["pert"])
            
    if ((yaml["frmType"] == "FORMATION_DEPUTY_ANCHOR") or
        (yaml["frmType"] == "FORMATION_DUAL_ABS")):
        # Define the deputy
        if yaml["deputy"]["stateType"] == "STATE_KEPEL":
            deputy = orb.Orbit(
                yaml["epoch"], 
                unpackKepEl(yaml["deputy"]["state"]),
                stateType = yaml["deputy"]["stateType"],
                pert = yaml["pert"]) 
            
    if ((yaml["frmType"] == "FORMATION_CHIEF_ANCHOR") or
        (yaml["frmType"] == "FORMATION_DEPUTY_ANCHOR")):
        # Define the relative state
        if yaml["relStateType"] == "RELSTATE_DOE":
            relState = unpackDoe(yaml["relState"]["state"])
        if ((yaml["relStateType"] == "RELSTATE_RECT_CLROE") or 
            (yaml["relStateType"] == "RELSTATE_CURV_CLROE")):
            relState = unpackClroe(yaml["relState"]["state"])
        if ((yaml["relStateType"] == "RELSTATE_RECTRIC") or 
            (yaml["relStateType"] == "RELSTATE_CURVRIC")):
            relState = unpackRelPosVel(yaml["relState"]["state"])
            
    return frm.Formation(
        chief, 
        deputy, 
        relState,
        yaml["frmType"],
        yaml["relStateType"],
        yaml["pert"])

def parseWptTbl(yaml, formation = None):
    """
    Parse Waypoint Plan yaml file

    Parameters
    ----------
    yaml : dictionary from yaml
    formation

    Returns
    -------
    wptTbl object

    """
    generator = yaml.pop("generator")
    if generator["type"] == "CLROE":
        wptTbl = wt.waypointTable()
        wptTbl.genClroeWpts(
            unpackClroe(generator["L"]),
            generator["meanMotion"],
            yaml["tblStartTime"],
            yaml["tblEndTime"],
            yaml["numWpts"])
        
        if yaml["relStateNatural"]:
            # Replace the formation deputy position with the natural path.
            # Assume curvilinear for long-range accuracy
            formation = frm.Formation(
                formation.chief, 
                formation.deputy, 
                unpackClroe(generator["L"]),
                "FORMATION_CHIEF_ANCHOR",
                "RELSTATE_CURV_CLROE",
                formation.chief.pert)
            
    return wptTbl, formation