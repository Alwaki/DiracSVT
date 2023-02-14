"""
Project:        Shell evolution of the dirac equation
                
Authors:        Daniel Karlsson & Alexander Kiessling
                (2021-2023)

Description:    Functionality for shooting method solver,
                this is the main looping point to call solving methods
                and handle results.
"""

from util import *
from numerical_solver import *
from boundary_conditions import *

def run_1Dsolve(state, Scenario):
    """ Run entire solver routine to generate binding energy solutions. This will call
    dirac solver until convergence/termination.

    Args:
        state (string): 
        Scenario (int, optional): defines parameter selection and potentials. Defaults to 1.

    Returns:
        [float]: solutions in the form of B and a0
    """
    data, p1, p2, p3, mn, mp = load_data("data.xlsx", "parameters.xlsx")
    
    if Scenario == 1:
        pfix = p1
        p = list(p1.values())
    if Scenario == 2:
        pfix = p2
        p = list(p2.values())
    if Scenario == 3:
        pfix = p3
        p = list(p3.values())
        
    pdict = {}
    state, t = state[0], state[1]
    pctrange = np.array(sorted(np.arange(-40, 41, 4), key = abs))/1000
    exp = data[(data.Name == state) & (data.t == t)].Experiment.values[0]
    
    samea = True
    if Scenario == 3:
        samea = False

    B, a0, rvals, FGvals  = solve_dirac(p, pfix, pdict, data, mn, mp, state = state, t = t, B0 = exp, Scenario = Scenario, samea = samea)    
    
    return B, a0, rvals, FGvals

def solve_dirac(p, pfix, pdict, data, mn, mp, state, t, l = False, k = False, N= False, Z = False, B0 = False, tensorV = 0.0, Scenario = False, samea = False):
    """ Runs an iteration of solving the dirac equation

    Args:
        p ([float]): list of parameter values
        pfix (): parameter object
        pdict (dictionary): dictionary of parameters
        data (struct): data object
        mn (float): neutron mass
        mp (float): proton mass
        state (string): state description
        t (bool): neutron or proton
        l (float, optional): quantum number. Defaults to False.
        k (int, optional): quantum number. Defaults to False.
        N (bool, optional): derived from atomic number. Defaults to False.
        Z (bool, optional): atomic number. Defaults to False.
        B0 (float/bool, optional): initial value B. Defaults to False.
        tensorV (float, optional): tensor parameter. Defaults to 0.
        Scenario (int, optional): scenario choice of potential. Defaults to False.
        samea (bool, optional): sets if same a should be used. Defaults to False.

    Returns:
        [float]: solutions in the form of B and a0
    """
    StateData = data[(data.Name == state) & (data.t == t)]
    for n in range(len(list(pdict.values()))):
        pdict[list(pdict.keys())[n]] = p[n]
    params = {**pdict, **pfix}
    V, kappa, lamda, r0, a, r0ls = params['V0'], params['kappa'], params['lambda'], params['r0'], params['a'], params['Rls']
    if samea:
        als = params['a']
    else:
        als = params['als']
    if Scenario == 3:
        kappa_so = params['k_so']
    N, Z, A, a0, xmin, xmax, xmatch, a0 = StateData.N.values[0], StateData.Z.values[0], StateData.N.values[0] + StateData.Z.values[0], 314000, StateData.xmin.values[0], StateData.xmax.values[0], StateData.xmatch.values[0], StateData.a0_d.values[0]
    if l == False and k == False:
        l, k, j = spectr(state.split()[1])
    B = B0 
    h = 0.0001

    # Determine appropriate potential
    if t == 1: #Proton state
        m = mp
        sigmaV0 = V*(1 + kappa*(N-Z)/A)
        if Scenario == 1:
            dV0 = -lamda*sigmaV0 
        elif Scenario == 2:
            dV0 = -lamda*V*(1 - kappa*(N-Z)/A) 
        elif Scenario == 3:
            dV0 = -lamda*V*(1 - kappa_so*(N-Z)/A) 
    elif t == -1: #Neutron state
        m = mn
        sigmaV0 = V*(1 - kappa*(N-Z)/A)
        if Scenario == 1:
            dV0 = -lamda*sigmaV0 
        elif Scenario == 2:
            dV0 = -lamda*V*(1 + kappa*(N-Z)/A) 
        elif Scenario == 3:
            dV0 = -lamda*V*(1 + kappa_so*(N-Z)/A)
    sigmaR = r0*(A)**(1/3)
    dR = r0ls*A**(1/3)
    sigmaa = a
    da = als
    
    error = 1
    iterations = 0
    while error > 0.0001 and iterations < 40:
        iterations += 1
        if B < 0:
            Foutbc, Goutbc, Finbc, Ginbc = BC(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m) #Set boundary conditions for B
        else:
            Foutbc, Goutbc, Finbc, Ginbc = BC_positive(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m) #Set boundary conditions for B positive
        inFG, rvals, FGvals = RK_integrate(xmax, xmatch, Finbc, Ginbc, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)
        outFG, rvalsout, FGvalsout = RK_integrate(xmin, xmatch, Foutbc, Goutbc, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)
        rvals = np.concatenate((rvals, rvalsout), axis = None)
        FGvals = np.concatenate((FGvals, FGvalsout), axis = 1)
        B1 = B*(1 + h)
        if B1 < 0:
            Foutbc, Goutbc, Finbc, Ginbc = BC(a0, xmin, xmax, l, B1, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m) #Set boundary conditions for B+h
        else:
            Foutbc, Goutbc, Finbc, Ginbc = BC_positive(a0, xmin, xmax, l, B1, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m) #Set boundary conditions for B positive
        dBinFG, drvals, dFGvals = RK_integrate(xmax, xmatch, Finbc, Ginbc, k, m, B1, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)
        dBoutFG, drvals, dFGvals = RK_integrate(xmin, xmatch, Foutbc, Goutbc, k, m, B1, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)

        DgfDB = (np.longdouble(Gap(dBoutFG, dBinFG))- np.longdouble(Gap(outFG, inFG)))/(B*h)
        da0outgf = outFG*(1+h) 
        try:
            DgfDa0 = (Gap(da0outgf, inFG) - Gap(outFG, inFG))/(a0*h)
        except FloatingPointError as e:
            print('DgfDa0 ERROR', e)
            return 100, 0, rvals, FGvals
        for i in range(len(DgfDa0)):
            for j in range(len(DgfDa0[i])):
                if DgfDa0[i][j] == 0:
                    DgfDa0[i][j] += 1e-12
        M = np.concatenate((DgfDB, DgfDa0)).reshape(2,2).T
        Old = np.array([B, a0])
        Cold = (np.dot(M, Old) - Gap(outFG,inFG)).reshape(2,1)
        try:
            New = np.linalg.solve(M, Cold)
            #New = np.linalg.lstsq(M, Cold, rcond = None)[0] #alternative way to solve, does not work well for most states
        except np.linalg.LinAlgError as e:
            print(e)
            print('Linalg error in main')
            return 100, 0, rvals, FGvals
        B = New[0][0]
        a0 = float(New[1][0])
        error = np.sqrt(np.linalg.norm(New - Old))
        if np.abs(B - Old[0]) < 1E-07:
            error = 0
        if B < -100 or B > 100:
            error = 0
            return 100, 0, rvals, FGvals
    return B, a0, rvals, FGvals