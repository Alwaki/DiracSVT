from util import *
from shooting_method import *

# Set state according to datasheet
state = ['208Pb 3p1/2', -1]

# Solver is run
B, a0 = run_1Dsolve(state, Scenario=1)
print("B: ", B, " a0: ", a0)