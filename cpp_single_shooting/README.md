## Table of contents
* [Setup](#setup)
* [Running](#running)
	
## Setup
The code requires compilation with C++17 or later, due to the dependency on bessel functions from the standard library. Moreover, there is also a dependecy on Eigen, which is available for download [here](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download).

## Running
Compilation has been tested using g++. The easiest way to compile and run is using vs code, for which settings have been included. Otherwise compilation can be performed through:

    compile instructions

after which the resulting binary can be executed through

    run instructions 

upon which the user will be prompted to select a scenario and a state. The corresponding states can be seen in the included data sheet.
These are read from two .csv files (data and parameters), which can be edited freely if desired to test with other parameters. An easy way to edit the files is to use the excel functionality to export to csv (excel versions are available in matlab and python directories).