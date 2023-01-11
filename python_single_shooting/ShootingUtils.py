import re
from numba import jit
import sys
import scipy.special as bessel
import numpy as np
from datetime import date
from datetime import datetime
import gvar as gv
from numpy.lib.scimath import sqrt as csqrt
import copy
import DiracData

#######################################################
#
#                 Physics functions
#
#######################################################

@jit(nopython = True)
def Tensor_fit(x, FG, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, isospin):
    """ Calculates the Dirac equation.
    :param x:
    :param FG:
    :param k:
    :param m:
    :param B:
    :param sigmaV0:
    :param dV0:
    :param sigmaR:
    :param dR:
    :param sigmaa:
    :param da:
    :param Z:
    :param tensorV:
    :param isospin:
    :returns:
    """
    dfgW = np.zeros((2,1))
    if isospin == 1:
        if x > sigmaR:
            Col = 0.0072923*Z/x
        else:
            Col = 0.0072923*Z*(3*sigmaR**2 - x**2)/(2*sigmaR**3)
    else:
        Col = 0
    sigma = sigmaV0/(1+np.exp((x-sigmaR)/sigmaa)) + Col                                     #V+S
    delta = dV0/(1+np.exp((x-dR)/da)) + Col                                                 #V-S
    U = tensorV/(1+np.exp((x-sigmaR)/sigmaa))                                               #tensor potential
    #Up = -(tensorV/sigmaa)*np.exp((x-sigmaR)/sigmaa)/(1+np.exp((x-sigmaR)/sigmaa))**2      #tensor potential derivative
    #dfgW[0] = (2*m+B-delta)*FG[1]+(U-k/x)*FG[0]
    #dfgW[1] = (-B + sigma)*FG[0]+(k/x-U)*FG[1]
    dfgW[1] = (2*m + B-delta)*FG[0]+(U-k/x)*FG[1]
    dfgW[0] = (-B + sigma)*FG[1]+(k/x-U)*FG[0]
    return dfgW

@jit(nopython = True)
def int_n_Tensor(xstart, xend, iniF, iniG, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, isospin):
    """ Attempts to solve the guessed initial value problem with RK4. 

    :param xstart:
    :param xend:
    :param iniF:
    :param iniG:
    :param k:
    :param m:
    :param B:
    :param sigmaV0:
    :param dV0:
    :param sigmaR:
    :param dR:
    :param sigmaa:
    :param da:
    :param Z:
    :param tensorV:
    :param isospin:
    :returns:
    """
    outputFG = np.zeros((2,1))
    h = 0.001
    rvals = np.zeros(len(range(1, int(1/h) + 1)))
    FGvals = np.zeros((2, len(range(1, int(1/h) + 1))))
       
    step = (xend - xstart)*h
    outputFG[0] = iniF
    outputFG[1] = iniG
    for i in range(1, int(1/h) + 1):
        x = (i-1)*step+xstart
        k1 = Tensor_fit(x, outputFG, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, isospin)
        k2 = Tensor_fit(x + step/2, outputFG + k1*step/2, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, isospin)
        k3 = Tensor_fit(x + step/2, outputFG + k2*step/2, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, isospin)
        k4 = Tensor_fit(x + step, outputFG + k3*step, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, isospin)
        outputFG += (k1 + 2*k2 + 2*k3 + k4)*step/6
        rvals[i - 1] = x
        FGvals[0][i - 1] = outputFG[0][0]
        FGvals[1][i - 1] = outputFG[1][0]
    return outputFG, rvals, FGvals

def spectr(state):
    """ Extracts the relevant quantum numbers of the given state

    :param state:   shell, example: 1p1/2 (string)
    :returns:       quantum numbers (float 3-tuple)
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

def BC(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m):
    if k < 0:
        Foutbc = -a0*xmin**(l+2)*(-B + sigmaV0/(1+np.exp(-sigmaR/sigmaa)))/(l+2-k)
        Goutbc = a0*xmin**(l+1)
    else:
        Foutbc = a0*xmin**l
        Goutbc = a0*xmin**(l+1)*(2*m+B-dV0/(1+np.exp(-dR/da)))/(l+k+1)    
    try:
        miu = np.sqrt(-2*m*B-B**2)
        Finbc = -np.sqrt(-B/(2*m + B))*np.exp(-miu*xmax)
    except FloatingPointError:
        miu = 1e-10
        Finbc = -np.sqrt(-B/(1e-10))*np.exp(-miu*xmax)
    Ginbc = np.exp(-miu*xmax)
    return Foutbc, Goutbc, Finbc, Ginbc

def BC_pos(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m, qm = 1):
    if k < 0:
        Foutbc = -a0*xmin**(l+2)*(-B + sigmaV0/(1+np.exp(-sigmaR/sigmaa)))/(l+2-k)
        Goutbc = a0*xmin**(l+1)
    else:
        Foutbc = a0*xmin**l
        Goutbc = a0*xmin**(l+1)*(2*m+B-dV0/(1+np.exp(-dR/da)))/(l+k+1)    
    if k < 0:
        Finbc, Ginbc = (bessel.spherical_jn(l, xmax) + bessel.spherical_yn(l, xmax)), np.sqrt(B**2 + 2*B*m)/(B + 2*m)*(bessel.spherical_jn(l + 1, xmax) + bessel.spherical_yn(l + 1, xmax)) #*spinor_sphereharm(l, k, j, qm, B) *spinor_sphereharm(l + 1, k, j, qm, B)
    else:
        Finbc, Ginbc = (bessel.spherical_jn(l, xmax) + bessel.spherical_yn(l, xmax)), np.sqrt(B**2 + 2*B*m)/(B + 2*m)*(bessel.spherical_jn(l - 1, xmax) + bessel.spherical_yn(l - 1, xmax)) #*spinor_sphereharm(l, k, j, qm, B) *spinor_sphereharm(l - 1, k, j, qm, B)
    norm = np.sqrt(np.sum(np.array([bessel.spherical_jn(l, xmax), bessel.spherical_yn(l, xmax)]))**2)
    return Foutbc, Goutbc, Finbc/norm, Ginbc/norm

def Gap(x, y):
    """ Calculates the difference between two numerical arrays

    :param x:       numpy array
    :param y:       numpy array
    :returns:       difference as numpy array
    """
    return np.array(x - y).T

def test_converge(B, Exp):
    """ Determines if method has diverged numerically.
        This is done by comparing the calculated binding energy with experimentally 
        known values, and seeing if the difference exceeds even a jump in shell levels.

    :param B:       numerically found binding energy value (float)
    :param Exp:     expected binding energy from experimental values (float)
    :returns:       convergence status (bool)
    """
    if Exp < 0:
        if B > Exp - 20 and B < 0:
            return True
    elif Exp > 0:
        if B < Exp + 20 and B > 0:
            return True
    return False

def run_1Dsolve(p, state, pfix, pdict, Scenario = 1):
    data = DiracData.data
    state, t = state[0], state[1]
    pctrange = np.array(sorted(np.arange(-40, 41, 4), key = abs))/1000
    exp = data[(data.Name == state) & (data.t == t)].Experiment.values[0]
    c = 0
    x = []
    y = []
    for b in range(0,20):
        B0 = exp + b*np.random.randn()/10
        B, a0 = solve_dirac(p, pfix, pdict, state = state, t = t, B0 = B0, plotWF = False, Scenario = Scenario, samea = True)    
        converge = test_converge(B, exp)
        if converge:
            break    
    paradjust = False
    for l in pctrange:
        if converge:
            break
        for i in range(len(p)):
            pars = copy.copy(p)
            #pars.iloc[i] = pars.iloc[i]*(1 + l)
            pars[i] = pars[i]*(1 + l/10)
            for b in range(0,20,4):
                B0 = exp + b*np.random.randn()/10
                B, a0 = solve_dirac(p, pfix, pdict, state = state, t = t, B0 = B0, plotWF = False, Scenario = Scenario, samea = False)    
                converge = test_converge(B, exp)
                if converge:
                   break    
            x.append(c)
            y.append(B)
            c+=1
            if B == -100:
                print('Not converged: ', l)
                print(list(pdict.keys())[i])
                #input('check')
                continue
            if converge:
                convpar = i
                parchange = l
                paradjust = True
                break
    if paradjust and parchange*p[convpar] != 0:
        print('Converged on parameter: ', list(pdict.keys())[convpar], ': Parameter Changed by: ', parchange*p[convpar])
        print('Old parameter: ', p[convpar], ' New parameter: ', pars[convpar])
    return B, a0

def solve_dirac(p, pfix, pdict, state, t, l = False, k = False, N= False, Z = False, B0 = False, tensorV = 0, plotWF = False, Scenario = False, samea = False, df = False ):
    
    # Load parameters
    mn = DiracData.mn
    mp = DiracData.mp
    StateData = DiracData.data[(DiracData.data.Name == state) & (DiracData.data.t == t)]
    for n in range(len(list(pdict.values()))):
        pdict[list(pdict.keys())[n]] = p[n]
    params = {**pdict, **pfix}
    V, kappa, lamda, r0, a, r0ls = params['V0'], params['kappa'], params['lambda'], params['r0'], params['a'], params['Rls']
    if samea:
        als = params['a']
    else:
        als = params['als']
    if Scenario == 3:
        kappa_so = params['kappa_so']
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
    while error > 0.0001:
        iterations += 1
        if B < 0:
            Foutbc, Goutbc, Finbc, Ginbc = BC(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m) #Set boundary conditions for B
        else:
            Foutbc, Goutbc, Finbc, Ginbc = BC_pos(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m, qm = 0) #Set boundary conditions for B positive
        inFG, rvals, FGvals = int_n_Tensor(xmax, xmatch, Finbc, Ginbc, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)
        outFG, rvalsout, FGvalsout = int_n_Tensor(xmin, xmatch, Foutbc, Goutbc, k, m, B, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)
        rvals = np.concatenate((rvals, rvalsout), axis = None)
        FGvals = np.concatenate((FGvals, FGvalsout), axis = 1)
        B1 = B*(1 + h)
        if B1 < 0:
            Foutbc, Goutbc, Finbc, Ginbc = BC(a0, xmin, xmax, l, B1, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m) #Set boundary conditions for B+h
        else:
            Foutbc, Goutbc, Finbc, Ginbc = BC_pos(a0, xmin, xmax, l, B1, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m, qm = 0) #Set boundary conditions for B positive
        dBinFG, drvals, dFGvals = int_n_Tensor(xmax, xmatch, Finbc, Ginbc, k, m, B1, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)
        dBoutFG, drvals, dFGvals = int_n_Tensor(xmin, xmatch, Foutbc, Goutbc, k, m, B1, sigmaV0, dV0, sigmaR, dR, sigmaa, da, Z, tensorV, t)

        DgfDB = (np.longdouble(Gap(dBoutFG, dBinFG))- np.longdouble(Gap(outFG, inFG)))/(B*h)
        da0outgf = outFG*(1+h) 
        try:
            DgfDa0 = (Gap(da0outgf, inFG) - Gap(outFG, inFG))/(a0*h)
        except FloatingPointError as e:
            print('DgfDa0 ERROR', e)
            return -100, 0
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
            return -100, 0
        B = New[0][0]
        a0 = float(New[1][0])
        error = np.sqrt(np.linalg.norm(New - Old))
        if np.abs(B - Old[0]) < 1E-07:
            error = 0
        if B < -100 or B > 100:
            #print('Not converging')
            #print(DgfDB)
            #input('stop')
            error = 0
            return -100, a0
        #else:
            #print('Converging')
            #print(DgfDB)
            #input('stop')
        #print(B)
    #print(rvals, FGvals)
    if plotWF:
        plotFG(df, rvals, FGvals, n, t)
    return B, a0

#######################################################
#
#                 Deprecated functions
#
#######################################################

def insert_error_value(B):
    value = 0 # Make Gvar centerered around B with uncertainty X, gaussisk fördelning 
    return value

def spinor_sphereharm(l, k, j, m, B, theta = 0, phi = 0):
    if k < 0:
        Y = 1/(2*l + 1)*np.array([np.sqrt(l + m + 1/2)*bessel.sph_harm(m - 1/2, l, theta, phi), np.sqrt(l - m + 1/2)])*bessel.sph_harm(m + 1/2, l, theta, phi)
    else:
        Y = 1/(2*l + 1)*np.array([-np.sqrt(l - m + 1/2)*bessel.sph_harm(m - 1/2, l, theta, phi), np.sqrt(l + m + 1/2)])*bessel.sph_harm(m + 1/2, l, theta, phi)
    return Y

def BC_2D(a0,c0, xmin, xmax, Omega, B, Sigma, sigmaR, sigmaa, dV0, dR, da, k, l, m):
    if k < 0:
        gm = a0*xmin**(l+2)
        fp = a0*((Omega + 1/2)/(B - Sigma))*xmin**(l+1)
        gp = c0*xmin**(l+2)
        fm = -c0*((Omega-1/2)/(B - Sigma))*xmin**(l+1)
    else:
        gm = a0*xmin**(l+2)
        fp = a0*((Omega + 1/2)/(B - Sigma))*xmin**(l+1)
        gp = c0*xmin**(l+2)
        fm = -c0*((Omega-1/2)/(B - Sigma))*xmin**(l+1)
    
    Ginbc = np.exp(-miu*xmax)
    Finbc = -np.sqrt(-B/(2*m + B))*np.exp(-miu*xmax)
    return Foutbc, Goutbc, Finbc, Ginbc

def run_solve(n, p, df, Scenario):
    pctrange = np.array(sorted(np.arange(-40, 41, 4), key = abs))/1000
    x = []
    y = []
    c = 0
    converge = False
    for b in range(0,20):
        B0 = df.iloc[n].Experiment + b*np.random.randn()/10
        B, a0 = solve_dirac(n, p, B0 = B0, plotWF = False, Scenario = Scenario, df = df, samea = True)
        converge = test_converge(B, df.iloc[n].Experiment)
        if converge:
            break    
    paradjust = False
    for l in pctrange:
        if converge:
            break
        for i in range(len(p)):
            pars = copy.copy(p)
            #pars.iloc[i] = pars.iloc[i]*(1 + l)
            pars[i] = pars[i]*(1 + l)
            for b in range(0,20,4):
                B, a0 = solve_dirac(n, pars, B0 = df.iloc[n].Experiment + b*np.random.randn()/10, plotWF = False, Scenario = Scenario, df = df, samea = True)
                converge = test_converge(B, df.iloc[n].Experiment)
                if converge:
                    #print(b, l)
                    break    
                #else:
                    #print(n)
            x.append(c)
            y.append(B)
            c+=1
            if B == -100:
                print('Not converged: ', l)
                continue
            if converge:
                convpar = i
                parchange = l
                paradjust = True
                break
    if not converge:
        error = df.iloc[n].Experiment + 20*np.random.randn()
        print('Inserted Error:', error)
        print(p)
        #input('press button')
        return error, 0, False, False, False
    #print('Binding Energy: ', B)
    #print('Experiment: ',df.iloc[n].Experiment)
    #print('dE : ', df.iloc[n].Experiment - B)
    if paradjust:
        print('Converged on parameter: ', convpar, ': Parameter Changed by: ', parchange*p[convpar])
        print('Old parameter: ', p[convpar], ' New parameter: ', pars[convpar])
        print('Converged on parameter: ', convpar, ': Parameter Changed by: ', parchange*p[convpar])
        print('Old parameter: ', p[convpar], ' New parameter: ', pars[convpar])
        return B , a0, convpar, parchange, pars[convpar]
    else:
        return B, a0, False, False, False
    #if len(x) > 1:
    #   plt.scatter(x,y)
    #   plt.show()    

#######################################################
#
#                 Book-keeping functions
#
#######################################################

def save_parameters(Pdf, plist, t, rms, filename, Scenario, samea): #TODO:Create dictionary instead to make it easier with change of fitparameters
    print("Used save_parameters")
    savelist = [str(Scenario)]
    tensorV = 0
    for x in plist:
        savelist.append(x)
    if Scenario == 3:
        if samea:
            savelist.insert(-1 ,0)
    while len(savelist) < len(Pdf.columns) - 2:
        savelist.append(0)
    savelist.extend([tensorV, rms])

    #Pdf.insert(8, 'k_so', '')
    Pdf.loc[datetime.now()] = savelist
    Pdf.to_excel(filename, header = True)
   

