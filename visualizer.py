# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 22:04:28 2025

@author: Hakim Lachnani
"""

import os

from plotting import plotDynamics as pltDyn
from plotting import plotNavigation as pltNav

def visAll(log, path, tag, settings):
    """
    Plot all available plots and animations

    """
    plotAll(log, path, tag, settings)
    animAll(log, path, tag, settings)
    
def plotAll(log, path, tag, settings):
    """
    Plot all available plots

    """
    pltDyn.plotAll(log, os.path.join(path, "dynamics"), tag, settings)
    if (settings["fsw"]["status"] == True):
        pltNav.plotAll(log, os.path.join(path, "navigation"), tag, settings)
    
def animAll(log, path, tag, settings):
    """
    Plot all available animations

    """
    pltDyn.animAll(log, os.path.join(path, "dynamics"), tag, settings)