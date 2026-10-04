# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import os

from plotting import plotDynamics as pltDyn

def visAll(log, path, tag, settings):
    """
    Plot all available plots and animations

    """
    pltDyn.plotAll(log, os.path.join(path, "dynamics"), tag, settings)
    pltDyn.animAll(log, os.path.join(path, "dynamics"), tag, settings)
    
def plotAll(log, path, tag, settings):
    """
    Plot all available plots

    """
    pltDyn.plotAll(log, os.path.join(path, "dynamics"), tag, settings)
    
def animAll(log, path, tag, settings):
    """
    Plot all available animations

    """
    pltDyn.animAll(log, os.path.join(path, "dynamics"), tag, settings)