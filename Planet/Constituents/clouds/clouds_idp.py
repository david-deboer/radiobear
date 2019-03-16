from __future__ import absolute_import, division, print_function
import math
import cmath
GHz = 29.9792458     # conversion from cm^-1 to GHz


def alpha(freq, T, P, cloud, cloud_dict, otherPar, units='dBperkm', path='./', verbose=False):
    """Adapted from Imke's code, but used Ulaby, Moore and Fung (see e.g. p310)."""
    alpha_cloud = []
    for f in freq:
        alpha = 0.0
        k = (2.0 * math.pi * f / GHz)  # *otherPar['refr']  # in cm^-1
        if 'ice' in otherPar.keys() and otherPar['ice'] > 0.0:
            fraction = cloud[cloud_dict['H2O']] / 0.9   # g/cm^3 / g/cm^3
            e = water(f, T)
            alpha += acloud(k, fraction, e)
        if 'water' in otherPar.keys() and otherPar['water'] > 0.0:
            fraction = cloud[cloud_dict['SOLN']] / 1.0
            e = water(f, T)
            alpha += acloud(k, fraction, e)
        if 'nh4sh' in otherPar.keys() and otherPar['nh4sh'] > 0.0:
            fraction = cloud[cloud_dict['NH4SH']] / 1.2
            nImke = 1.7 - 0.05j
            nDeBoer = 1.74 - 0.001j
            n = 1.7 - 0.005j
            e = n**2
            alpha += acloud(k, fraction, e)
        if 'nh3ice' in otherPar.keys() and otherPar['nh3ice'] > 0.0:
            fraction = cloud[cloud_dict['NH3']] / 1.6
            nImke = 1.3 - 0.05j
            nDeBoer = 1.3 - 0.005j
            n = 1.3 - 0.0001j
            e = n**2
            alpha += acloud(k, fraction, e)
        if 'h2sice' in otherPar.keys() and otherPar['h2sice'] > 0.0:
            fraction = cloud[cloud_dict['H2S']] / 1.5
            nImke = 1.3 - 0.01j
            nDeBoer = 1.15 - 0.001j
            n = 1.15 - 0.0001j
            e = n**2
            alpha += acloud(k, fraction, e)
        if 'ch4' in otherPar.keys() and otherPar['ch4'] > 0.0:
            fraction = cloud[cloud_dict['CH4']] / 1.0
            n = 1.3 - 0.00001j
            e = n**2
            alpha += acloud(k, fraction, e)
        if alpha < 0.0:
            print("Warning:  cloud alpha<0.  Reset to 0")
            alpha = 0.0
        if units == 'dBperkm':
            alpha *= 434294.5
        alpha_cloud.append(alpha)
    if verbose:
        print('cloud absorption: ', alpha_cloud)
    return alpha_cloud


def acloud(k, fraction, e):
    """See Ulaby, Moore and Fung p310 or Janssen et al 102"""
    K = (e - 1.0) / (e + 2.0)
    alp = 3.0 * k * fraction * (-K.imag)
    return alp


FR = [1.0E8, 3.0E8, 1.0E9, 2.0E9, 3.0E9, 5.0E9, 1.0E10, 3.0E10, 1.0E11]
EIMAG = [8.0E-3, 1.5E-3, 8.0E-4, 1.0E-3, 1.2E-3, 1.5E-3, 3.0E-3, 8.0E-3, 2.0E-2]


def water(freq, T):
    T_Celsius = T - 273.0
    freqHz = freq * 1.0E9
    if T_Celsius >= 0.0:
        RelT = 1.1109E-10 - T_Celsius * 3.824E-12 + (T_Celsius**2) * 6.938E-14 - (T_Celsius**3) * 5.096E-16
        E0 = 88.045 - 0.4147 * T_Celsius + (T_Celsius**2) * 6.295E-4 + (T_Celsius**3) * 1.075E-5
        if E0 < 0.0:
            E0 = 0.0
        EINF = 4.9
        E1 = EINF + (E0 - EINF) / (1.0 + (freqHz * RelT)**2)
        E2 = freqHz * RelT * (E0 - EINF) / (1.0 + (freqHz * RelT)**2)
        if E2 < 0.0:
            E2 = 0.0
    else:
        E1 = 3.17
        LF = math.log10(freqHz)
        for j in range(len(FR) - 1):
            if FR[j + 1] >= freqHz:
                break
        LF0 = math.log10(FR[j])
        LF1 = math.log10(FR[j + 1])
        DLF = (LF - LF0) / (LF1 - LF0)
        X0 = math.log10(EIMAG[j])
        X1 = math.log10(EIMAG[j + 1])
        DX = X0 + DLF * (X1 - X0)
        E2 = 10.0**DX
    dielectricConstant = E1 - E2 * 1.0j
    refractiveIndex = cmath.sqrt(dielectricConstant)
    return dielectricConstant
