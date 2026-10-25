# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:09:44 2025

@author: Hakim Lachnani
"""

import numpy as np
import simulator

def scenario(sim, dvTbl, acc = None):
    """
    Runs the chiefDv scenario, which executes a (series of) chief Delta-V(s)

    Parameters
    ----------
    sim : Simulator object
        see simulator.py
    dvTbl : delta-v table object
        see dvTable.py
    acc : acceleration capability in km/s^2

    Returns
    -------
    None.

    """
    print("SCR: Running chiefDv for {} seconds".format(sim.settings["simDuration"]))
    
    # Target the subsequent waypoints
    for ii in range(dvTbl.t):
        if dvTbl.t[ii] < sim.settings["simDuration"]:
            # Propagate to burn epoch
            sim.run(dvTbl.t[ii])
            # Compute burn in Eci
            dvEci = dvTbl.asInr(0, sim.frm.dcmInr2Ric)
            # Execute burn
            if acc == None:
                sim.impulsiveBurn(dvEci, sv = "Chief", enablePert = False)
            else:
                # Scale the burn and acceleration to meet the DV
                dvAcc = acc
                accEciDir = dvEci/np.linalg.norm(dvEci)
                sim.frm.deputy.aCtrlInEci = accEciDir * dvAcc
    
    sim.run(sim.settings["simDuration"])
    sim.terminate()
    


