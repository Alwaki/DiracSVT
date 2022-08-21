
import numpy as np
import sys
import time
import timeit
np.set_printoptions(threshold = sys.maxsize)
import pandas as pd
np.set_printoptions(precision=3, threshold=sys.maxsize)
from ShootingUtils import *
import DiracData
import scipy.optimize as op

def getsign(p):
    if p < 0:
        return 1
    else:
        return -1


state = ['16O 1p1/2', -1]

data = DiracData.data
p1 = DiracData.p1
p = {'V0': p1['V0'], 'r0': p1['r0'], 'Rls': p1['Rls']}
pnames = list(p1.keys())
pfixed = {pnames[pnames.index(x)]:p1[x] for x in pnames if x not in list(p.keys())}

B, a0 = run_1Dsolve(list(p1.values()), state, {}, p1, Scenario = 1)
print(B, a0)
    
