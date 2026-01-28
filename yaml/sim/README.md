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

The dynamics dt controls the integration step, and therefore accuracy of the truth solution.

The formation parameters should be nominally enabled, but may be disabled to improve simulation speed.