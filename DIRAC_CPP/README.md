## Table of contents
* [Setup](#setup)
* [Running](#running)
* [Notes](#notes)
	
## Setup
The code requires compilation with C++17 or later, due to the dependency on special functions from the standard library. Moreover, there is also a dependecy on Eigen, which is available for download [here](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download). If using linux (ubuntu), eigen can be installed using the following steps

    sudo apt update && sudo apt upgrade

    sudo apt install libeigen3-dev

## Running
A finished binary is included, which can be run on linux using 

    ./main

Compilation has been tested using g++ (GNU). The easiest way to compile and run is using vs code. If running vs code, simply add the C/C++ extension, and configure it for c++17 (or later). Otherwise compilation can be performed through:

    g++ --std c++17 -g main.cpp -o main

after which the resulting binary can be executed as above.

## Notes
- When running the binary, the user will be prompted to select a scenario and a state. These are read from two .csv files (data and parameters), which can be edited freely if desired to test with other parameters. An easy way to edit the files is to use the excel functionality to export to csv (excel versions are available in matlab and python directories).
- The C++ version does not support plotting the resulting wavefunction, but instead exports these to a .txt document. This document is formatted as three columns, with each column using a comma and a space as delimiting characters. The columns are radial distance, f-component and g-component, respectively.