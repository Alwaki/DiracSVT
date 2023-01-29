## Table of contents
* [Setup](#setup)
* [Running](#running)

## Setup
This version makes use of several python libraries. These can be installed as a batch using pip, following

    pip install scipy pandas numba numpy

Note that it might be necessary to use pip3, depending on your python environment. Moreover other package handlers are also usable, such as anaconda.

## Running

Generally, the only file one should need to interact with to run the code is the main.py file. In the file one can select which state to investigate (see data sheet for which ones are available), and then run by

    python main.py

or alternatively modify the main file to an executable and run it directly. Note that if executing through python, it might be necessary to call the file through python3 instead.