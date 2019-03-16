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
    filename = os.path.join(path, 'h2s.npz')
    if verbose:
        print("Reading h2s lines from {}".format(filename))
    global data
    data = np.load(filename)


def alpha(freq, T, P, X, P_dict, otherPar, units='dBperkm', path='./', verbose=False):
    """Computes absorption due to h2s"""

    # Read in data if needed
    if data is None:
        readInputFiles(path, verbose)

    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_h2s = P * X[P_dict['H2S']]
    f0 = data['f0']
    I0 = data['I0']
    E = data['E']
    GH2S = data['GH2S']
    GH2 = 1.960
    GHe = 1.200
    ZH2 = 0.000
    ZHe = 0.000
    ZH2S = 0.000
    C = 1.0
    D = 1.28
    n_dvl = 0.7
    n_int = 3.0 / 2.0
    delta = D * P_h2s
    gamma = np.power((T0 / T), n_dvl) * (GH2 * P_h2 + GHe * P_he + GH2S * P_h2s)
    g2 = gamma**2
    zeta = gamma
    z2 = zeta**2
    ITG = I0 * np.exp(-((1.0 / T) - (1.0 / T0)) * E * hck)

    alpha_h2s = []
    for f in freq:
        f2 = f**2
        num = (gamma - zeta) * f2 + (gamma + zeta) * (np.power(f0 + delta, 2.0) + g2 - z2)
        den = np.power((f2 - np.power(f0 + delta, 2.0) - g2 + z2), 2.0) + 4.0 * f2 * g2
        shape = GHz * 2.0 * np.power(f / f0, 2.0) * num / (np.pi * den)
        alpha_h2s.append(np.sum(shape * ITG))

    alpha_h2s = coef * (P_h2s / T0) * pow((T0 / T), n_int + 2) * np.array(alpha_h2s)
    if units == 'dBperkm':
        alpha_h2s *= 434294.5

    return alpha_h2s
