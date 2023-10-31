## Table of contents
* [Setup](#setup)
* [Running](#running)


## Setup
This version makes use of several python libraries. These can be installed as a batch using pip, following

    pip install scipy matplotlib pandas numba numpy==1.20.3 openpyxl

If you already have a specific version which is incompatible, you can force an install of a specific version as follows:
    
    pip install --force-reinstall numpy==1.20.3

## Running

Generally, the only file one should need to interact with to run the code is the main.py file. In the file one can select which state to investigate (see data sheet for which ones are available), and then run by executing the file from terminal from the scripts directory. An example of passing arguments can be seen below

    python3 main.py --state 7 --plot False --scenario 1 

State 7 in this case corresponds to 40Ca 1f7/2 (see the data tabular sheet for other states). The expected output from this example would be approximately B: -9.48. If the solver becomes severely numerically unstable for a case, it will return B:  100  a0:  0.

If it is desired to run all available states at once, then one can pass 0 as an argument for state.

