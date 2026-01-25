# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:38:52 2025

@author: Hakim Lachnani
"""

import os

from tkinter import filedialog
import pickle
        
def logLoad(path = None):
    """
    Load DATB Log File

    Returns
    -------
    None.

    """
    if path is None:
        path = filedialog.askopenfilename(defaultextension = '.pkl',
                                               initialdir = os.getcwd())
    
    with open(path, 'rb') as inp:
        return pickle.load(inp)
    