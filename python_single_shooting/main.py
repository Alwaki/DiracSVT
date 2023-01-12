#!/usr/bin/env python3

"""
Project:        Shell evolution of the dirac equation with 
                tensor potential
                
Authors:        Alexander W. Kiessling & Daniel Karlsson 
                (2021-2023)

Description:    Main executable file of program
"""

from util import *
from shooting_method import *

# Set state according to datasheet
state = ['132Sn 2f7/2', -1]

# Solver is run
B, a0 = run_1Dsolve(state, Scenario=1)
print("B: ", B, " a0: ", a0)