#!/usr/bin/env python

"""
setup.py file for SWIG dynamics utilities
"""

from distutils.core import setup, Extension
import numpy

dynamicsUtils_module = Extension('_dynamicsUtils',
                           sources=['dynamicsUtils_wrap.c', 'dynamicsUtils.c'],
                           include_dirs =  [numpy.get_include(),'.']
                           )

setup (name = 'dynamicsUtils',
       version = '0.1',
       author      = "Hakim Lachnani",
       description = """DATB Dynamics Utilities in C""",
       ext_modules = [dynamicsUtils_module],
       py_modules = ["dynamicsUtils"],
       )
