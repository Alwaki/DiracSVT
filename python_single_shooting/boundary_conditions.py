"""
Project:        Shell evolution of the dirac equation
                
Authors:        Alexander W. Kiessling & Daniel Karlsson 
                (2021-2023)

Description:    Boundary conditions for the radial dirac 
                equation.

"""

import scipy.special as bessel
import numpy as np

def BC(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m):
    """ Calculates positive boundary condition

    Args:
        a0 (float): diffusivity of nuclear surface
        xmin (float): min range position for potential
        xmax (float): max range position for potential
        l (float): quantum number
        B (float): energy
        sigmaV0 (_type_): sigmaV0 (float): potential well depth descriptor
        sigmaR (float): potential well depth descriptor
        sigmaa (float): potential parameter
        dV0 (float): potential well depth descriptor
        dR (float): range of potential well
        da (float): potential parameter
        k (float): quantum number
        m (float): mass

    Returns:
        [float]: boundary conditions at separate ends for f and g
    """
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

def BC_positive(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m):
    """ Calculates positive boundary condition

    Args:
        a0 (float): diffusivity of nuclear surface
        xmin (float): min range position for potential
        xmax (float): max range position for potential
        l (float): quantum number
        B (float): energy
        sigmaV0 (_type_): sigmaV0 (float): potential well depth descriptor
        sigmaR (float): potential well depth descriptor
        sigmaa (float): potential parameter
        dV0 (float): potential well depth descriptor
        dR (float): range of potential well
        da (float): potential parameter
        k (float): quantum number
        j (float): quantum number
        m (float): mass

    Returns:
        [float]: boundary conditions at separate ends for f and g
    """
    if k < 0:
        Foutbc = -a0*xmin**(l+2)*(-B + sigmaV0/(1+np.exp(-sigmaR/sigmaa)))/(l+2-k)
        Goutbc = a0*xmin**(l+1)
    else:
        Foutbc = a0*xmin**l
        Goutbc = a0*xmin**(l+1)*(2*m+B-dV0/(1+np.exp(-dR/da)))/(l+k+1)    
    if k < 0:
        Finbc, Ginbc = (bessel.spherical_jn(l, xmax) + bessel.spherical_yn(l, xmax)), np.sqrt(B**2 + 2*B*m)/(B + 2*m)*(bessel.spherical_jn(l + 1, xmax) + bessel.spherical_yn(l + 1, xmax)) 
    else:
        Finbc, Ginbc = (bessel.spherical_jn(l, xmax) + bessel.spherical_yn(l, xmax)), np.sqrt(B**2 + 2*B*m)/(B + 2*m)*(bessel.spherical_jn(l - 1, xmax) + bessel.spherical_yn(l - 1, xmax))
    norm = np.sqrt(np.sum(np.array([bessel.spherical_jn(l, xmax), bessel.spherical_yn(l, xmax)]))**2)
    return Foutbc, Goutbc, Finbc/norm, Ginbc/norm