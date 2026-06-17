#!/usr/bin/env python

"""
setup.py file for SWIG kinematics utilities
"""

from distutils.core import setup, Extension
import numpy

kinematicsUtils_module = Extension('_kinematicsUtils',
                           sources=['kinematicsUtils_wrap.c', 'kinematicsUtils.c'],
                           include_dirs =  [numpy.get_include(),'.']
                           )

setup (name = 'kinematicsUtils',
       version = '0.1',
       author      = "Hakim Lachnani",
       description = """DATB Kinematics Utilities in C""",
       ext_modules = [kinematicsUtils_module],
       py_modules = ["kinematicsUtils"],
       )
