from util import *
from shooting_method import *

state = ['208Pb 3p1/2', -1]
stuff = load_data("data.xlsx", "parameters.xlsx")
data = stuff[0]
for i in range(0, 5):
    state = [data['Name'][i], -1]
    B, a0 = run_1Dsolve(state, Scenario=1)
    print(B)
    #print(a0)
    #print("B: ", B, " a0: ", a0)