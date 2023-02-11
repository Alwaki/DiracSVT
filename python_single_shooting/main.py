#!/usr/bin/env python3

"""
Project:        Shell evolution of the dirac equation
                
Authors:        Daniel Karlsson & Alexander Kiessling
                (2021-2023)

Description:    Main executable file of program
"""

from util import *
from shooting_method import *

# Set state according to datasheet
#state = ['16O 1p1/2', -1]
d = pd.read_excel(os.getcwd()+"\\python_single_shooting\\"+"data.xlsx")
s = d["Name"].values
t = d["t"].values
for i, x in enumerate(t):
    
    state = [s[i], x]

    # Solver is run
    B, a0, rvals, FGvals = run_1Dsolve(state, Scenario=2)
    print(B)


#this is a change
# Print out results and plot wavefunction
#print("B: ", B, " a0: ", a0)
#plotWF(rvals, FGvals, state)