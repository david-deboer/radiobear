from __future__ import absolute_import, division, print_function
import os.path
import numpy as np

# Some constants
coef = 7.244E+21     # coefficient from GEISA_PP.TEX eq. 14
T0 = 296.0           # reference temperature in K
hck = 1.438396       # hc/k  [K cm]
GHz = 29.9792458     # conversion from cm^-1 to GHz

data = None


def readInputFiles(path, verbose):
    filename = os.path.join(path, 'co.npz')
    if verbose:
        print("Reading co lines from {}".format(filename))
    global data
    data = np.load(filename)


def alpha(freq, T, P, X, P_dict, otherPar, units='dBperkm', path='./', verbose=False):
    # ##Voigt coefficients from Janssen p67
    PLimits = [0.001, 0.1]
    avoigt = [122.60793178, 214.38238869, 181.92853309, 93.15558046, 30.18014220, 5.91262621, 0.56418958, 0.0]
    bvoigt = [122.60793178, 352.73062511, 457.33447878, 348.70391772, 170.35400182, 53.99290691, 10.47985711, 1.0]

    # Read in data if needed
    if data is None:
        readInputFiles(path, verbose)

    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_co = P * X[P_dict['CO']]
    f0 = data['f0']
    I0 = data['I0']
    E = data['E']
    GH2 = 1.960
    GHe = 1.200
    GCO = 6.000
    n_dvl = 0.7
    n_int = 3.0 / 2.0
    gamma = pow((T0 / T), n_dvl) * (GH2 * P_h2 + GHe * P_he + GCO * P_co)
    g2 = gamma**2
    zeta = 0.0
    z2 = zeta**2
    delta = 0.0
    ITG = I0 * np.exp(-((1.0 / T) - (1.0 / T0)) * E * hck)
    w = (P - PLimits[0]) / (PLimits[1] - PLimits[0])
    if w < 0.0:
        w = 0.0
    elif w > 1.0:
        w = 1.0

    alpha_co = []
    for f in freq:
        f2 = f**2
        shape_Voigt = np.zeros(len(f0))
        if P <= PLimits[1] or otherPar == 'voigt' or otherPar == 'diff':
            # ##Doppler broadening Janssen p59
            betaD = 4.3e-7 * np.sqrt(T / 28.0) * f
            # ##Voigt Janssen p67
            num = np.zeros(len(f0), dtype='complex128')
            den = np.zeros(len(f0), dtype='complex128')
            xi = gamma / betaD + (1.0j) * (f - f0) / betaD
            for jjj in range(len(avoigt)):
                num += avoigt[jjj] * (xi**jjj)
                den += bvoigt[jjj] * (xi**jjj)
            val = num / den
            shape_Voigt = GHz * (1.0 / (np.sqrt(np.pi) * betaD)) * val.real
        shape_VVW = np.zeros(len(f0))
        if P >= PLimits[0] or otherPar == 'vvw' or otherPar == 'diff':
            num = (gamma - zeta) * f2 + (gamma + zeta) * (np.power(f0 + delta, 2.0) + g2 - z2)
            den = np.power((f2 - np.power(f0 + delta, 2.0) - g2 + z2), 2.0) + 4.0 * f2 * g2
            shape_VVW = GHz * 2.0 * np.power(f / f0, 2.0) * num / (np.pi * den)
        shape = w * shape_VVW + (1.0 - w) * shape_Voigt
        if otherPar == 'voigt':
            alpha_co.append(np.sum(shape_Voigt))
        elif otherPar == 'vvw':
            alpha_co.append(np.sum(shape_VVW))
        elif otherPar == 'diff':
            alpha_co.append(np.sum(shape_Voigt - shape_VVW))
        else:
            alpha_co.append(np.sum(shape * ITG))

    alpha_co = coef * (P_co / T0) * pow((T0 / T), n_int + 2) * np.array(alpha_co)
    if units == 'dBperkm':
        alpha_co *= 434294.5

    return alpha_co
