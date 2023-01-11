import re, os, sys, copy
import pandas as pd
import numpy as np
import math
import scipy.constants as constants

def spectr(state):
    """ Extracts the relevant quantum numbers of the given state

    Args:
        state (string): shell, example: 1p1/2

    Returns:
        [float]: quantum numbers
    """
    lnot = [('s', 0), ('p', 1), ('d', 2), ('f', 3), ('g', 4), ('h', 5), ('i', 6), ('j', 7), ('k', 8)]
    l = int([ls[1] for ls in lnot if ls[0] == ''.join(re.findall('\\D',state.replace('/','')))][0])
    j = float(state[2:].split('/')[0])/float(state[2:].split('/')[1])
    if j == l + 1/2:
        k = -(l + 1)
    elif j == l - 1/2:
        k = l 
    else:
        print("Error jkl")
        sys.exit()
    return l, k, j

def load_data(filename1, filename2):
    """ Loads the necessary project data

    Args:
        file1 (string): excel file for data
        file2 (string): excel file for parameters

    Returns:
        list: parameter data, scenario specific data, mass
    """
    data = pd.read_excel(os.getcwd()+"\\python_single_shooting\\"+filename1)
    pardf= pd.read_excel(os.getcwd()+"\\python_single_shooting\\"+filename2)
    data.index = data['Unnamed: 0']
    data.pop('Unnamed: 0')
    p1 = {pardf.columns[x]:pardf[pardf.Scenario == 1].values[0][x] for x in range(1,8)}
    p2 = {pardf.columns[x]:pardf[pardf.Scenario == 2].values[0][x] for x in range(1,7)}
    p3 = {pardf.columns[x]:pardf[pardf.Scenario == 3].values[0][x] for x in range(1,8)}

    meq = "mass energy equivalent in MeV"
    mn = constants.physical_constants["neutron " + meq][0] #energy in MeV
    mp = constants.physical_constants["proton " + meq][0] #energy in MeV
    return data, p1, p2, p3, mn, mp

def Gap(x, y):
    """ Calculates the difference between two numerical arrays

    Args:
        x (array): A
        y (array): B

    Returns:
        array: difference between A and B
    """
    return np.array(x - y).T

def test_converge(B, Exp):
    """ Determines if method has diverged numerically.
    This is done by comparing the calculated binding energy with experimentally 
    known values, and seeing if the difference exceeds even a jump in shell levels.

    Args:
        B (float): numerically found binding energy value
        Exp (float): expected binding energy from experimental values

    Returns:
        bool: convergence status
    """
    return math.isclose(B, Exp, rel_tol=0.1)