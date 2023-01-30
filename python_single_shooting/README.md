## Table of contents
* [Setup](#setup)
* [Running](#running)
* [Notes](#notes)


## Setup
This version makes use of several python libraries. These can be installed as a batch using pip, following

    pip3 install scipy pandas numba numpy openpyxl

Note that it might be necessary to use pip, depending on your python environment (same for python calls). Moreover other package handlers are also usable, such as anaconda.

## Running

Generally, the only file one should need to interact with to run the code is the main.py file. In the file one can select which state to investigate (see data sheet for which ones are available), and then run by

    python3 main.py

or alternatively modify the main file to an executable and run it directly.

## Notes
- If the solver diverges, it will return values of B: -100, a0: 0. 
- The longdouble in shooting_method.py is deprecated in later editions of numpy. Changing this to a float64 will allow the code to run, but will reduce numerical stability and produce slightly different results in cases with tendencies towards divergence.
- The current file pathing is different between linux and windows, see the load_data function in util.py. For windows it works to use double backslashes, on linux this will work with single forward slashes.
