# Init
DATB initialization YAML files define the initial state and perturbations of the formation and are structured as follows:
- epoch: J2000 epoch time in seconds
- chief: chief orbit (may be omitted by setting to 0)
    - state: chief orbit state
    - stateType: chief orbit state type
- deputy: deputy orbit (may be omitted by setting to 0)
    - state: deputy orbit state
    - stateType: deputy orbit state type
- relState: (may be omitted by setting to 0)
    - state: relative state
    - relStateType: relative state type
- frmType: formation type
- pert: perturbation settings
    - jnum: Number of J-gravitational parameters
    - solarGrav: enable/disable flag for Solar gravity
    - lunarGrav: enable/disable flag for lunar gravity
    - SRP: enable/disable flag for solar radiation pressure
    - drag: enable/disable flag for atmospheric drag

Do define the formation, two of the following are required:
- chief state
- deputy state
- relative state

frmType will indicate which two define the formation.