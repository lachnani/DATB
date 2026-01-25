# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 18:33:02 2025

@author: Hakim Lachnani
"""

import numpy as np
from dynamics import formation as frm
from dynamics import orbit as orb

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
            chiefState = np.array((
                yaml["chief"]["state"]["a"],
                yaml["chief"]["state"]["e"],
                yaml["chief"]["state"]["i"]*orb.D2R,
                yaml["chief"]["state"]["RAAN"]*orb.D2R,
                yaml["chief"]["state"]["argP"]*orb.D2R,
                yaml["chief"]["state"]["M"]*orb.D2R))
            chief = orb.Orbit(
                yaml["epoch"], 
                chiefState,
                stateType = yaml["chief"]["stateType"],
                pert = yaml["pert"])
            
    if ((yaml["frmType"] == "FORMATION_DEPUTY_ANCHOR") or
        (yaml["frmType"] == "FORMATION_DUAL_ABS")):
        # Define the deputy
        if yaml["deputy"]["stateType"] == "STATE_KEPEL":
            deputyState = np.array((
                yaml["deputy"]["state"]["a"],
                yaml["deputy"]["state"]["e"],
                yaml["deputy"]["state"]["i"]*orb.D2R,
                yaml["deputy"]["state"]["RAAN"]*orb.D2R,
                yaml["deputy"]["state"]["argP"]*orb.D2R,
                yaml["deputy"]["state"]["M"]*orb.D2R))
            deputy = orb.Orbit(
                yaml["epoch"], 
                deputyState,
                stateType = yaml["deputy"]["stateType"],
                pert = yaml["pert"]) 
            
    if ((yaml["frmType"] == "FORMATION_CHIEF_ANCHOR") or
        (yaml["frmType"] == "FORMATION_DEPUTY_ANCHOR")):
        # Define the relative state
        if yaml["relStateType"] == "RELSTATE_DOE":
            relState = np.array((
                yaml["relState"]["state"]["da"],
                yaml["relState"]["state"]["de"],
                yaml["relState"]["state"]["di"]*orb.D2R,
                yaml["relState"]["state"]["dRAAN"]*orb.D2R,
                yaml["relState"]["state"]["darg"]*orb.D2R,
                yaml["relState"]["state"]["dM"]*orb.D2R))
        if ((yaml["relStateType"] == "RELSTATE_RECT_CLROE") or 
            (yaml["relStateType"] == "RELSTATE_CURV_CLROE")):
            relState = np.array((
                yaml["relState"]["state"]["A0"],
                yaml["relState"]["state"]["alpha"]*orb.D2R,
                yaml["relState"]["state"]["xOff"],
                yaml["relState"]["state"]["yOff"],
                yaml["relState"]["state"]["B0"],
                yaml["relState"]["state"]["beta"]*orb.D2R))
        if ((yaml["relStateType"] == "RELSTATE_RECTRIC") or 
            (yaml["relStateType"] == "RELSTATE_CURVRIC")):
            relState = np.array((
                yaml["relState"]["state"]["relPosR"],
                yaml["relState"]["state"]["relPosI"],
                yaml["relState"]["state"]["relPosC"],
                yaml["relState"]["state"]["relVelR"],
                yaml["relState"]["state"]["relVelI"],
                yaml["relState"]["state"]["relVelC"]))
            
    return frm.Formation(
        chief, 
        deputy, 
        relState,
        yaml["frmType"],
        yaml["relStateType"],
        yaml["pert"])