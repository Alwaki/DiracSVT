"""
Project:        Shell evolution of the dirac equation 
                
Authors:        Daniel Karlsson & Alexander Kiessling
                (2021-2023)

Description:    Numerical methods to solve initial value problem through
                integration. Methods are decorated with jit handle,
                which reduces debug value. Remove jit decorator
                in case of debug.
"""

from numba import jit
import numpy as np

@jit(nopython = True)
def point_solver(x, FG, k, m, B, sigmaV0, dV0, sigmaR, dR,
                 sigmaa, da, Z, tensorV, isospin):
    """ Calculation of slope and potentials at point x

    Args:
        x (float): position in potential
        FG ([float]): function values f and g
        k (float): quantum number related to l and j
        m (float): mass
        B (float): energy
        sigmaV0 (float): potential well depth descriptor
        dV0 (float): potential well depth descriptor
        sigmaR (float): potential well depth descriptor
        dR (float): range of potential well
        sigmaa (float): potential parameter
        da (float): potential parameter
        Z (float): atomic number
        tensorV (float): tensor potential parameter
        isospin (bool): neutron or proton 

    Returns:
        [float]: derivative values of f' and g'
    """
    dfgW = np.zeros((2,1))
    if isospin == 1:
        if x > sigmaR:
            Col = 0.0072923*Z/x
        else:
            Col = 0.0072923*Z*(3*sigmaR**2 - x**2)/(2*sigmaR**3)
    else:
        Col = 0
    sigma = sigmaV0/(1+np.exp((x-sigmaR)/sigmaa)) + Col                 # vector potential + scalar potential
    delta = dV0/(1+np.exp((x-dR)/da)) + Col                             # vector potential - scalar potential
    U = tensorV/(1+np.exp((x-sigmaR)/sigmaa))                           # tensor potential
    dfgW[1] = (2*m + B-delta)*FG[0]+(U-k/x)*FG[1]
    dfgW[0] = (-B + sigma)*FG[1]+(k/x-U)*FG[0]
    return dfgW

@jit(nopython = True)
def RK_integrate(xstart, xend, iniF, iniG, k, m, B, sigmaV0, dV0,
                 sigmaR, dR, sigmaa, da, Z, tensorV, isospin):
    """ Runge-Kutta integration

    Args:
        xstart (float): beginning position of integration
        xend (float): end position of integration
        iniF (float): initial value of f
        iniG (float): initial value of g
        k (float): quantum number related to l and j
        m (float): mass
        B (float): energy
        sigmaV0 (float): potential well depth descriptor
        dV0 (float): potential well depth descriptor
        sigmaR (float): potential well depth descriptor
        dR (float): range of potential well
        sigmaa (float): potential parameter
        da (float): potential parameter
        Z (float): atomic number
        tensorV (float): tensor potential parameter
        isospin (bool): neutron or proton 

    Returns:
        [float]: f and g values
    """
    outputFG = np.zeros((2,1))
    h = 0.001
    
    # Storage, not used atm
    rvals = np.zeros(len(range(1, int(1/h) + 1)))
    FGvals = np.zeros((2, len(range(1, int(1/h) + 1))))
       
    step = (xend - xstart)*h
    outputFG[0] = iniF
    outputFG[1] = iniG
    for i in range(1, int(1/h) + 1):
        x = (i-1)*step+xstart
        k1 = point_solver(x, outputFG, k, m, B, sigmaV0, dV0, sigmaR,
                          dR, sigmaa, da, Z, tensorV, isospin)
        k2 = point_solver(x + step/2, outputFG + k1*step/2, k, m, B,
                          sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z,
                          tensorV, isospin)
        k3 = point_solver(x + step/2, outputFG + k2*step/2, k, m, B,
                          sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z,
                          tensorV, isospin)
        k4 = point_solver(x + step, outputFG + k3*step, k, m, B,
                          sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z,
                          tensorV, isospin)
        outputFG += (k1 + 2*k2 + 2*k3 + k4)*step/6
        rvals[i - 1] = x
        FGvals[0][i - 1] = outputFG[0][0]
        FGvals[1][i - 1] = outputFG[1][0]
    return outputFG, rvals, FGvals