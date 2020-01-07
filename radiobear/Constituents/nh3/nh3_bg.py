import math
import os.path
import numpy as np
from radiobear.Constituents import parameters

# Some constants
coef = 7.244E+21     # coefficient from GEISA_PP.TEX eq. 14
T0 = 296.0           # reference temperature in K
hck = 1.438396       # hc/k  [K cm]
GHz = 29.9792458     # conversion from cm^-1 to GHz


data = None


def readInputFiles(par):
    """This reads in the data files for nh3"""
    filename = os.path.join(par.path, 'nh3.npz')
    if par.verbose:
        print("Reading nh3 lines from {}".format(filename))
    global data
    data = np.load(filename)


def alpha(freq, T, P, X, P_dict, other_dict, **kwargs):

    par = parameters.setpar(kwargs)
    # Read in data if needed
    global data
    if data is None:
        readInputFiles(par)

    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_nh3 = P * X[P_dict['NH3']]
    f0 = data['f0']
    I0 = data['I0']
    E = data['E']
    G0 = data['G0']
    nlin = len(f0)

    GH2 = 2.318
    GHe = 0.790
    GNH3 = 0.750
    ZH2 = 1.920
    ZHe = 0.300
    ZNH3 = 0.490
    # C = 1.0075 + (0.0308 + 0.552 * P_h2 / T) * P_h2 / T
    D = -0.45
    n_dvl = 2.0 / 3.0
    n_int = 3.0 / 2.0
    delta = D * P_nh3

    alpha_nh3 = []
    for f in freq:
        f2 = f**2
        alpha = 0.0
        for i in range(nlin):
            gamma = pow((T0 / T), n_dvl) * (GH2 * P_h2 + GHe * P_he + G0[i] * GNH3 * P_nh3)
            g2 = gamma**2
            zeta = pow((T0 / T), n_dvl) * (ZH2 * P_h2 + ZHe * P_he + G0[i] * ZNH3 * P_nh3)
            z2 = zeta**2
            ITG = I0[i] * math.exp(-((1.0 / T) - (1.0 / T0)) * E[i] * hck)
            num = (gamma - zeta) * f2 + (gamma + zeta) * (pow(f0[i] + delta, 2.0) + g2 - z2)
            den = pow((f2 - pow(f0[i] + delta, 2.0) - g2 + z2), 2.0) + 4.0 * f2 * g2
            shape = GHz * 2.0 * pow(f / f0[i], 2.0) * num / (math.pi * den)
            alpha += shape * ITG

        a = coef * (P_nh3 / T0) * pow((T0 / T), n_int + 2) * alpha
        if par.units == 'dBperkm':
            a *= 434294.5
        alpha_nh3.append(a)

    return alpha_nh3
