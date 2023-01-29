#!/usr/bin/env python3

"""
Project:        Shell evolution of the dirac equation
                
Authors:        Alexander W. Kiessling & Daniel Karlsson 
                (2021-2023)

Description:    Main executable file of program
"""

from util import *
from shooting_method import *

# Set state according to datasheet
state = ['16O 1p1/2', -1]

# Solver is run
B, a0, rvals, FGvals = run_1Dsolve(state, Scenario=1)

# Print out results and plot wavefunction
print("B: ", B, " a0: ", a0)
plotWF(rvals, FGvals, state)