## Table of contents
* [Setup](#setup)
* [Running](#running)
* [Notes](#notes)
	
## Setup
The code depends on the library Eigen, which is available for download [here](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download). If using linux (ubuntu), eigen can be installed using the following steps

    sudo apt update && sudo apt upgrade

    sudo apt install libeigen3-dev

## Running
A finished binary is included in the lib directory, which (provided user is in /lib directory) can be run on linux using 

    ./DIRAC_SOLVER

The user is then prompted to select a scenario and a state. These are in accordance with the accompanying data and parameter files, which can be inspected in the lib directory. The code will return the value of the binding energy B, and an auxiliary parameter a0. If the solver diverges, these will be returned as 100 and 0, respectively.

If a full build is desired, cmake files are provided. Assuming user is in the base directory of DIRAC_CPP, then Cmake can be invoked through

    cmake -S . -B lib/

This will produce make files in the lib directory. Navigating to /lib, we can then use

    make

which will produce a binary executable that can be run as shown above.

## Notes
- The C++ version does not support plotting the resulting wavefunction, but instead exports these to a .txt document. This document is formatted as three columns, with each column using a comma and a space as delimiting characters. The columns are radial distance, f-component and g-component, respectively.