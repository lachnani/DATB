# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 11:09:44 2025

@author: Hakim Lachnani
"""

from mathematics import mathematicsUtils as uMath

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
    
    # Execute DV table
    for ii in range(dvTbl.numDv):
        if dvTbl.t[ii] < sim.settings["simDuration"]:
            # Propagate to burn epoch
            sim.runUntil(dvTbl.t[ii])
            # Compute burn in Eci
            dvEci = dvTbl.asInr(ii, sim.frm.dcmInr2Ric)
            # Execute burn
            print(f"SCR: Executing [{1000*dvTbl.dv[0,ii]:.4f},{1000*dvTbl.dv[1,ii]:.4f},{1000*dvTbl.dv[2,ii]:.4f}] m/s {dvTbl.frame[ii]} DV")
            if acc == None:
                sim.impulsiveBurn(dvEci, sv = "Chief", enablePert = False)
            else:
                # Scale the burn and acceleration to meet the DV
                dvEciDir, dvMag = uMath.unitVec(dvEci)
                dvDur = dvMag/acc
                dvCyc = round(dvDur/sim.settings["dynamics"]["dt"])
                if dvCyc == 0:
                    # Acceleration is too large, execute impulsively
                    sim.impulsiveBurn(dvEci, sv = "Chief", enablePert = False)
                else:
                    dvDur = dvCyc * sim.settings["dynamics"]["dt"]
                    dvAcc = dvMag/dvDur
                    sim.frm.deputy.aCtrlInEci = dvEciDir * dvAcc
                    sim.runFor(dvDur)
                    sim.frm.deputy.aCtrlInEci = uMath.zeroVec()
    
    sim.runUntil(sim.settings["simDuration"])
    sim.terminate()
    


