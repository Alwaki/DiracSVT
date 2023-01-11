import scipy.special as bessel
import numpy as np

def BC(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, m):
    """_summary_

    Args:
        a0 (_type_): _description_
        xmin (_type_): _description_
        xmax (_type_): _description_
        l (_type_): _description_
        B (_type_): _description_
        sigmaV0 (_type_): _description_
        sigmaR (_type_): _description_
        sigmaa (_type_): _description_
        dV0 (_type_): _description_
        dR (_type_): _description_
        da (_type_): _description_
        k (_type_): _description_
        m (_type_): _description_

    Returns:
        _type_: _description_
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

def BC_positive(a0, xmin, xmax, l, B, sigmaV0, sigmaR, sigmaa, dV0, dR, da, k, j, m, qm = 1):
    """_summary_

    Args:
        a0 (_type_): _description_
        xmin (_type_): _description_
        xmax (_type_): _description_
        l (_type_): _description_
        B (_type_): _description_
        sigmaV0 (_type_): _description_
        sigmaR (_type_): _description_
        sigmaa (_type_): _description_
        dV0 (_type_): _description_
        dR (_type_): _description_
        da (_type_): _description_
        k (_type_): _description_
        j (_type_): _description_
        m (_type_): _description_
        qm (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: _description_
    """
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