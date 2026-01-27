# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:09:44 2025

@author: Hakim Lachnani
"""

import numpy as np
from planning import wptTbl as wt
from guidance import wptGuidance as wg
from guidance import linearizedModels as lm

def scenario(sim, relPosRic = None, rptFreq = None, dtStart = None, stm = None):
    """
    Runs the rptWpt scenario, which repeats a waypoint indefinitely.

    Parameters
    ----------
    sim : Simulator object
        see simulator.py
    relPosRic : repeat waypoint relative position
    rptFreq : repeat waypoint frequency (per period)

    Returns
    -------
    None.

    """
    print("SCR: Running rptWpt for {} seconds".format(sim.settings["simDuration"]))
    
    if relPosRic == None:
        # Prompt user for position
        print("SCR: User input repeated waypoint:")
        relPosRic = np.zeros((3,))
        relPosRic[0] = float(input("SCR: - Enter Radial component in km:"))
        relPosRic[1] = float(input("SCR: - Enter In-Track component in km:"))
        relPosRic[2] = float(input("SCR: - Enter Cross-Track component in km:"))
        
    if rptFreq == None:
        # Prompt user for frequency
        print("SCR: User input frequency:")
        rptFreq = float(input("SCR: - Enter per-period frequency:"))
        
    if dtStart == None:
        # Promp user for table start time
        print("SCR: User input table start time:")
        dtStart = float(input("SCR: - Enter start time in seconds:"))
    
    print("SCR: Creating waypoint table:")
    print(f"SCR: - Pos = [{relPosRic[0]:.2f},{relPosRic[1]:.2f},{relPosRic[2]:.2f}] km")
    print("SCR: - Freq = {} wpts per period".format(rptFreq))
    print("SCR: - Start Time = {} seconds".format(dtStart))
    wptTbl = wt.waypointTable()
    period = 2*np.pi/sim.frm.chief.meanMotion
    dt = period/rptFreq
    numWpts = int(np.ceil(sim.settings["simDuration"]/dt))
    wptTbl.genRptWpts(relPosRic, sim.t+dtStart, dt, numWpts)
    
    if stm == None:
        # Assume HCW wpt targeting
        stm = lm.hcw(sim.frm.chief.meanMotion)
    
    # Immediately burn to target the first waypoint
    dvRic = wg.wptBurn(stm.Phi(sim.frm.t,wptTbl.t[0]), 
                    sim.frm.relPosRectRic, 
                    sim.frm.relVelRectRic, 
                    wptTbl.relPosRic[:,0])
    dvEci = np.matmul(np.transpose(sim.frm.dcmInr2Ric), dvRic)
    sim.impulsiveBurn(dvEci)
    print(f"SCR: Executing [{1000*dvRic[0]:.4f},{1000*dvRic[1]:.4f},{1000*dvRic[2]:.4f}] m/s RIC DV to target waypoint 1")
    
    for wpt in range(numWpts-1):
        sim.run(wptTbl.t[wpt])
        dvRic = wg.wptBurn(stm.Phi(sim.frm.t,wptTbl.t[wpt+1]), 
                        sim.frm.relPosRectRic, 
                        sim.frm.relVelRectRic, 
                        wptTbl.relPosRic[:,wpt+1])
        dvEci = np.matmul(np.transpose(sim.frm.dcmInr2Ric), dvRic)
        sim.impulsiveBurn(dvEci)
        print(f"SCR: Executing [{1000*dvRic[0]:.4f},{1000*dvRic[1]:.4f},{1000*dvRic[2]:.4f}] m/s RIC DV to target waypoint {wpt+2}")
    
    sim.run(sim.settings["simDuration"])
    sim.terminate()
    


