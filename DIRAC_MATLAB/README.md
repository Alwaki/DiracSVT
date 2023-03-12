## Table of contents
* [Setup](#setup)
* [Running](#running)
	
## Setup
The matlab version does not require any additional toolboxes to run. As such it can be run out of the box without any setup.

## Running
Navigate to the main.m file, which is available in the DIRAC_MATLAB base directory. Running this script will perform the shooting method on the given state. The corresponding states can be seen in the included data sheet. If it is desired to run all states, set the state as 0. If the solver severely diverges for a case, it will return B = 100, a0 = 0.