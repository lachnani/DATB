# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:09:44 2025

@author: Hakim Lachnani
"""

def scenario(sim):
    """
    Runs the Idle scenario

    Parameters
    ----------
    sim : Simulator object
        see simulator.py

    Returns
    -------
    None.

    """
    print("SCR: Running idle for {} seconds".format(sim.settings["simDuration"]))
    
    sim.run(sim.settings["simDuration"])
    
    sim.terminate()
    
    
