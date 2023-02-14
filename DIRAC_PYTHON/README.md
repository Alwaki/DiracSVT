## Table of contents
* [Setup](#setup)
* [Running](#running)


## Setup
This version makes use of several python libraries. These can be installed as a batch using pip, following

    pip3 install scipy pandas numba numpy==1.20.3 openpyxl

If you already have a specific version which is incompatible, you can force an install of a specific version as follows:
    
    pip3 install --force-reinstall numpy==1.20.3

## Running

Generally, the only file one should need to interact with to run the code is the main.py file. In the file one can select which state to investigate (see data sheet for which ones are available), and then run by executing the file from terminal from the scripts directory. An example of passing arguments can be seen below

    python main.py --state "40Ca 1f7/2" --plot False --scenario 2 --particle -1

The expected output from this example would be B:  -9.76521264734975  a0:  9257096.937828392. If the solver becomes severely numerically unstable, it will return B:  100  a0:  0.

Note that the arguments regarding selection of particle and state need to match in accordance with the datasheet.

