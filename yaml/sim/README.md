# Sim
DATB simulation YAML files define simulation settings and are structured as follows:
- name: file name
- dynamics:
    - status: enable/disable flag
    - dt: time step for propagation
    - kov: keep-out-volume
- simDuration: simulation duration in seconds
- log: log settings
    - status: enable/disable flag
    - dt: time step for logging
- visualizer: visualizer settings
    - plots: enable/disable result plots
    - animations: enable/disable result animations
- fsw: FSW settings
    - status: enable/disable flag
    - dt: time step for FSW
- formation: formation settings
    - relStates: enable/disable flag for formation relative state computations
    - environments: enable/disable flag for formation environment computations
    - measurements: enable/disable flag for formation measurement computations
    - orbit: orbit settings
        - environments: enable/disable flag for orbit environement computations
        - elements: enable/disable flag for orbit elements computations

The dynamics dt controls the integration step, and therefore accuracy of the truth solution. Integration step size should be kept smaller than the following, per the guidelines from Shuster [[1]](#1):
- LEO: 1s
- GEO: 15s
- Molniya: 1s
- GTO: 1s

The formation parameters should be nominally enabled, but may be disabled to improve simulation speed.

## References
<a id="1">[1]</a> 
S. P. Shuster,
*A Survey and Performance Analysis of Orbit Propagators for LEO, GEO, and Highly Elliptical Orbits.*
phdthesis, Utah State University, May 2017.