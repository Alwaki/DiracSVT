## Table of contents
* [Setup](#setup)
* [Running](#running)
	
## Setup
The matlab version does not require any additional toolboxes to run. As such it can be run out of the box without any setup.

## Running
Navigate to the main.m file, which is available in the DIRAC_MATLAB base directory. Running this script will perform the shooting method on the given state. To change the evaluated state and scenario, the parameters 

	scenario = 1;
	element_state_index = 1;

can be changed to other values. The corresponding states can be seen in the included data sheet. If it is desired to run all states, set the state as 0. If the solver diverges for a case, it will return B = 100, a0 = 0.
