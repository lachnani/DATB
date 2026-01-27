# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:09:44 2025

@author: Hakim Lachnani
"""

import numpy as np

def scenario(sim):
    """
    Runs the halfTimeBurn scenario, which performs a 1 m/s RSO in-track burn 
    half way through the sim duration. Assumes 0.01 m/s^2 of acceleration such 
    that the burn takes 100s

    Parameters
    ----------
    sim : Simulator object
        see simulator.py

    Returns
    -------
    None.

    """
    print("SCR: Running halfTimeBurn for {} seconds".format(sim.settings["simDuration"]))
    
    idleTime = sim.settings["simDuration"]/2 - 50
    
    print("SCR: Idling for {} seconds".format(idleTime))
    sim.run(idleTime)
    
    print("SCR: Performing 1 m/s RSO In-Track Burn with 0.01 m/s^2")
    accEci = np.matmul(np.transpose(sim.frm.dcmInr2Ric), np.array([0,0.00001,0]))
    sim.frm.chief.aCtrlInEci = accEci # In km/s^2
    sim.run(idleTime + 100)
    print("SCR: RSO Burn Complete")
    
    print("SCR: Idling for {} seconds".format(idleTime))
    sim.frm.chief.aCtrlInEci = np.zeros((3,))
    sim.run(sim.settings["simDuration"])
    
    sim.terminate()
    
    
