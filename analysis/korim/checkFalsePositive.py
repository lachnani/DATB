# -*- coding: utf-8 -*-
"""
Created on Fri Aug  1 22:17:50 2025

@author: Hakim Lachnani
"""

import matplotlib.pyplot as plt
from itertools import product, combinations
import numpy as np

import simulator
import parser

# For the cases where KORIM fails, we want to show that the true trajectory 
# barely clips the edges of the KOV
print("MC: Checking False-Positive results ...")
opsaKorimFP = np.argwhere(opsaKorimErr > 0) 
unique, unique_indices = np.unique(opsaKorimFP[:,0], return_index = True)  
opsaKorimFP = opsaKorimFP[unique_indices,:]
numFP = np.shape(opsaKorimFP)[0]
opsaKorimFPMinRange = np.zeros((numFP,))
settings_fp = settings_case.copy()
settings_fp['log']['dt'] = 30
for othNum in range(numOth):
    oth = oths[othNum] 
    # Initialize plot
    korimFig = plt.figure()
    # plt.title(f'KORIM False-Positive Trajectories: {oth:.0f} Orbit OTH')
    ax = korimFig.add_subplot(projection="3d")
    origin = ax.scatter(0, 0, 0, color='b', linewidths = 3)
    ax.set_xlabel("Curvilinear In-Track [km]")
    ax.set_ylabel("Curvilinear Radial [km]")
    ax.set_zlabel("Curvilinear Cross-Track [km]")
    # Create cube
    r = [-kovLim[0], kovLim[0]]
    i = [-kovLim[1], kovLim[1]]
    c = [-kovLim[2], kovLim[2]]
    for s, e in combinations(np.array(list(product(i, r, c))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0] or np.sum(np.abs(s-e)) == i[1]-i[0] or np.sum(np.abs(s-e)) == c[1]-c[0]:
            ax.plot3D(*zip(s, e), color="b")
    for i in range(numFP):
        if othNum == opsaKorimFP[i,1]:
            simNum = opsaKorimFP[i,0]
            settings_fp["simDuration"] = oths[opsaKorimFP[i,1]]*orbPeriod
            # Rerun each false positive simulation with a 30s log rate
            print("MC: Re-running Sim {}".format(simNum+1))
            # Compute initial CLROEs using uniform distribution
            frm_sim = frm_base.copy()
            frm_sim["relState"]["state"]["A0"] = initConds[simNum,0]
            frm_sim["relState"]["state"]["alpha"] = initConds[simNum,1]
            frm_sim["relState"]["state"]["xOff"] = initConds[simNum,2]
            frm_sim["relState"]["state"]["yOff"] = initConds[simNum,3]
            frm_sim["relState"]["state"]["B0"] = initConds[simNum,4]
            frm_sim["relState"]["state"]["beta"] = initConds[simNum,5]
            
            # Initialize the sim
            sim = simulator.Simulator(settings_fp, parser.parseFormation(frm_sim), 0, quiet = True)
            
            # Run the Sim
            sim.run(simDuration)
            
            # Plot the result
            iPlot = sim.log.relPosCurvRic.loc['I',0:sim.log.i].to_numpy()
            iBreachInd = np.abs(iPlot) <= kovLim[1]
            rPlot = sim.log.relPosCurvRic.loc['R',0:sim.log.i].to_numpy()
            rBreachInd = np.abs(rPlot) <= kovLim[0]
            cPlot = sim.log.relPosCurvRic.loc['C',0:sim.log.i].to_numpy()
            cBreachInd = np.abs(cPlot) <= kovLim[2]
            
            allBreachInd = np.bitwise_and(np.bitwise_and(iBreachInd, rBreachInd), cBreachInd)
            
            if any(allBreachInd):
                ax.plot(iPlot[allBreachInd], 
                        rPlot[allBreachInd], 
                        cPlot[allBreachInd],
                        linewidth=1.0, color = 'r')[0]
                
                # Save the minimum range
                opsaKorimFPMinRange[i] = np.min(np.linalg.norm(np.c_[rPlot[allBreachInd],iPlot[allBreachInd],cPlot[allBreachInd]], axis=1))
            else:
                opsaKorimFPMinRange[i] = np.max(kovLim)

    # Clean up the plot
    plt.grid()
    plt.gca().invert_xaxis()
    ax.set_aspect('equal', adjustable = 'datalim')