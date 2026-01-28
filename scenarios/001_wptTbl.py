# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:09:44 2025

@author: Hakim Lachnani
"""

import numpy as np
import simulator
from guidance import wptGuidance as wg
from guidance import linearizedModels as lm

def scenario(sim, wptTbl, accVeh = None, stm = None, 
             burnComputeTimeScale = 1.5):
    """
    Runs the wptTbl scenario, which executes a waypoint table

    Parameters
    ----------
    sim : Simulator object
        see simulator.py
    wptTbl : waypoiny table object
        see wptTable.py
    accVeh : vehicle acceleration capability in km/s^2

    Returns
    -------
    None.

    """
    print("SCR: Running wptTbl for {} seconds".format(sim.settings["simDuration"]))
    
    if stm == None:
        # Assume HCW wpt targeting
        stm = lm.hcw(sim.frm.chief.meanMotion)
        
    # Compute the nominal burn table
    burnTblT, burnTblDv = wg.wptBurnTable(
        stm, wptTbl, 0, sim.frm.relPosRectRic, sim.frm.relVelRectRic)
    
    # Immediately burn to target the first waypoint
    dvEci = np.matmul(np.transpose(sim.frm.dcmInr2Ric), burnTblDv[:,0])
    sim.impulsiveBurn(dvEci)
    print(f"SCR: Executing [{1000*burnTblDv[0,0]:.4f},{1000*burnTblDv[2,0]:.4f},{1000*burnTblDv[2,0]:.4f}] m/s RIC DV to target waypoint 1")
    
    # Target the subsequent waypoints
    for wpt in range(wptTbl.numWpts):
        if accVeh == None:
            # Execute Implulsive burn
            sim.run(wptTbl.t[wpt])
            dvRic = wg.wptBurn(stm.Phi(sim.frm.t,wptTbl.t[wpt+1]), 
                            sim.frm.relPosRectRic, 
                            sim.frm.relVelRectRic, 
                            wptTbl.relPosRic[:,wpt+1])
            dvEci = np.matmul(np.transpose(sim.frm.dcmInr2Ric), dvRic)
            sim.impulsiveBurn(dvEci)
        else:
            # Estimate the burn time from the nominal burn table
            dvRicNom = burnTblDv[:,wpt+1]
            dvTimeNom = np.linalg.norm(dvRicNom)/accVeh
            # Simulate to burn decision time
            sim.run(wptTbl.t[wpt] - dvTimeNom*burnComputeTimeScale/2)
            # Propagate current states to waypoint time
            rWpt, vWpt = stm.prop(sim.frm.t, 
                            sim.frm.relPosRectRic, 
                            sim.frm.relVelRectRic,
                            wptTbl.t[wpt])
            dvRic = wg.wptBurn(stm.Phi(wptTbl.t[wpt],wptTbl.t[wpt+1]), 
                            rWpt, 
                            vWpt, 
                            wptTbl.relPosRic[:,wpt+1])
            dvRic = simulator.pertVec(dvRic, sim.burnPertAngDeg, sim.burnPertMagPerc)
            # Re-estimate the burn time
            dvMag = np.linalg.norm(dvRic)
            dvTime = dvMag/accVeh
            # Simulate to burn start time
            sim.run(wptTbl.t[wpt] - dvTime/2)
            # Execute Burn
            dvCyc = np.ceil(dvTime/sim.settings["dynamics"]["dt"])
            dvAcc = dvMag/dvCyc*sim.settings["dynamics"]["dt"]
            if dvAcc > accVeh:
                print("SCR: Error in Delta-V Quantization!")
            accRicDir = dvRic/dvMag
            for ii in range(dvCyc):
                accEciDir = np.matmul(np.transpose(sim.frm.dcmInr2Ric), accRicDir)
                sim.frm.deputy.aCtrlInEci = accEciDir * dvAcc
                sim.run(sim.settings["dynamics"]["dt"])
            sim.frm.deputy.aCtrlInEci = np.zeros((3,))
        print(f"SCR: Executing [{1000*dvRic[0]:.4f},{1000*dvRic[1]:.4f},{1000*dvRic[2]:.4f}] m/s RIC DV to target waypoint {wpt+2}")
    
    sim.run(sim.settings["simDuration"])
    sim.terminate()
    


