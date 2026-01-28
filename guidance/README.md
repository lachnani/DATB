# Guidance
The DATB Guidance module contains models and guidance tools to determine desired relative states.

## Linearized Models
Linearized models of relative motion are contained in linearizedModels.py. Models include:
- HCW (Hill-Clohessy-Wiltshire)
- CLROE (Circular Linearized Relative Orbit Elements)
- YA (Yamanaka-Ankersen)
- GA (Gim-Alfriend)
- sGA (Simplified Gim-Alfriend)
- doe (differential orbit elements)
- dee (differential equinoctial elements)
- dAmico (D'Amico relative orbit elements)
- dAmico2 (D'Amico relative orbit elements plus drag)
- ELROE (Elliptical Linearized Relative Orbit Elements)

Each model is an independent class which may have unique attributes. All model classes possess the following methods:
- A: state matrix A from the dynamics equation xDot = A*x + B*u
- B: input matrix B from the dynamics equation xDot = A*x + B*u
- Phi: state transition matrix 
- prop: propagates state

## Waypoint Guidance
Tools for guidance between position defined waypoints are contained in wptGuidance.py. Functions include:
- twoImpulseBurn: RIC frame burn vectors for two-impulse rendezvous
- wptBurn: RIC frame burn vector for targeting a position waypoint
- wptBurnTable: computes a series of waypoint burns for a waypoint table





