"""
Title:        		    DiracSVT
                
Authors:        	    Alexander Kiessling, Daniel Karlsson, 
			            Yuxin Zhao, Chong Qi

Version:		        1.0 (03/2023)	

Project Description:    Numerical solution of the Dirac equation with scalar,
			            vector and tensor potentials

File Description:	    Utility functions for data loading and 
                        less clearly defined tools.

"""

import re, os, sys, copy, argparse
import pandas as pd
import numpy as np
import math
import scipy.constants as constants
import matplotlib.pyplot as plt

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
    try:
        # Windows pathing
        cur_path = os.path.dirname(__file__)
        new_path = os.path.relpath('..\\data\\', cur_path)
        data = pd.read_excel(new_path+"\\"+filename1)
        pardf= pd.read_excel(new_path+"\\"+filename2)
    except:
        # Linux pathing
        cur_path = os.path.dirname(__file__)
        new_path = os.path.relpath('../data/', cur_path)
        data = pd.read_excel(new_path+"/"+filename1)
        pardf= pd.read_excel(new_path+"/"+filename2)
        
    data.index = data['Unnamed: 0']
    data.pop('Unnamed: 0')
    p1 = {pardf.columns[x]:pardf[pardf.Scenario == 1].values[0][x] for x in range(1,9)}
    p2 = {pardf.columns[x]:pardf[pardf.Scenario == 2].values[0][x] for x in range(1,8)}
    p3 = {pardf.columns[x]:pardf[pardf.Scenario == 3].values[0][x] for x in range(1,9)}

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

def plotWF(rvals, FGvals, state):
    """ Plots wavefunction

    Args:
        rvals ([float]): radial component
        FGvals ([float]): split wavefunction (two component solution)
        state ([tuple]): tuple containing string and integer. String is used as title.
    """
    x, y = zip(*sorted(zip(rvals,FGvals[0])))
    plt.plot(x, y, c = 'r', label = 'f')
    x, y = zip(*sorted(zip(rvals,FGvals[1])))
    plt.plot(x, y, c = 'b', label = 'g')
    plt.legend()
    plt.title(state[0])
    plt.xlabel('r')
    plt.ylabel('Wavefunction')
    plt.show()