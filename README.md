# DATB: Differential Astrodynamics Test Bench

DATB (Differential Astrodynamics Test Bench) is a modular planning, test, simulation, and visualization suite for relative orbital motion.

DATB is built with a C backend and a Python front end for easy scripting and rapid runtimes.

The main components of DATB are:
- Dynamics for orbit and formation propagation.
- Planning for trajectory design and optimization.
- Guidance for trajectory models and propagation algorithms.
- Estimation for simulation of onboard filter algorithms.
- Control for simulation of onboard controls algorithms.

The analysis folder contains scripts and functions for specialized analyses, including Monte Carlo simulations.

The yaml folder contains consifuguration scripts.

The scenarios folder contains mission and test scenarios.

Created by Hakim Lachnani.

## Dependencies
- yaml
- tkinter
- SWIG
- ffmpeg
