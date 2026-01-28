# Planning
The DATB planning module contains tools for planning RPO missions. The include waypoint table generation and optimization, as well as passive safety tools.

## Waypoint Table
The waypoint table class in contained in wptTbl.py. Waypoint table attributes are:
- numWpts: number of waypoints in the table
- t: time of waypoints in the table
- relPosRic: RIC frame relative position of the waypoints

The basic methods of the class are appendWpts() for adding waypoints to the table and increment() for incrementing the table and removing the first (now passed) waypoint. There are also several methods which simplify waypoint generation:
- genClroeWpts: generates waypoints along natural trajectory defined by CLROEs
- genRptWpts: generates repeated waypoints
- genFmcWpts: generates forced motion circumnavigation (FMC) waypoints

## Passive Safety
passiveSafety.py contains functions and tools for evaluating passive safety. These include:
- checkKovBreach: Determines whether the relative position is within the keep-out-volume (KOV)
- clroeProp: operational passive safety assessment (OPSA) using CLROE propagation
- stmProp: OPSA using state transition matrix (STM) propagation from a linearized model
- korim: OPSA using keep-out-region intersection method (KORIM)