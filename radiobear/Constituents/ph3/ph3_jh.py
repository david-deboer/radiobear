from __future__ import absolute_import, division, print_function
import os.path
import numpy as np

# Some constants
coef = 7.244E+21        # coefficient from GEISA_PP.TEX eq. 14
T0 = 300.0              # reference temperature in K
hck = 1.438396          # hc/k  [K cm]
GHz = 29.9792458        # conversion from cm^-1 to GHz
PI = np.pi

# Set data arrays and values
data = None
data_wgt = {}


def readInputFiles(path, verbose=False):
    filename = os.path.join(path, 'ph3jh.npz')
    if verbose:
        print("Reading ph3 lines from  {}".format(filename))
    global data
    data = np.load(filename)
    filename = os.path.join(path, 'PH3WGT.npz')
    if verbose:
        print("Reading ph3 wgts from {}".format(filename))
    global data_wgt
    data_wgt = np.load(filename)


def alpha(freq, T, P, X, P_dict, otherPar, units='dBperkm', path='./', verbose=False):
    """Computes the absorption due to ph3"""

    # Read in data if needed
    if data is None:
        readInputFiles(path, verbose=verbose)

    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_ph3 = P * X[P_dict['PH3']]
    f0 = data['f0']
    I0 = data['I0']
    E = data['E']
    WgtI0 = data_wgt['WgtI0']
    WgtFGB = data_wgt['WgtFGB']
    WgtSB = data_wgt['WgtSB']

    # Line parameters
    GH2 = 3.2930
    GHe = 1.6803
    GPH3 = 4.2157
    n_dvl = 2.0 / 3.0
    n_int = 3.0 / 2.0
    zeta = 0.0
    z2 = zeta**2
    delta = 0.0
    gamma = np.power((T0 / T), n_dvl) * (GH2 * P_h2 + GHe * P_he) * WgtFGB + np.power((T0 / T), 1.0) * GPH3 * P_ph3 * WgtSB
    g2 = gamma**2
    ITG = I0 * WgtI0 * np.exp(-((1.0 / T) - (1.0 / T0)) * E * hck)

    alpha_ph3 = []
    for f in freq:
        f2 = f**2
        num = (gamma - zeta) * f2 + (gamma + zeta) * ((f0 + delta)**2 + g2 - z2)
        den = (f2 - (f0 + delta)**2 - g2 + z2)**2 + 4.0 * f2 * g2
        shape = GHz * 2.0 * ((f / f0)**2) * num / (PI * den)
        alpha_ph3.append(np.sum(shape * ITG))

    alpha_ph3 = coef * (P_ph3 / T0) * np.power((T0 / T), n_int + 2) * np.array(alpha_ph3)
    if units == 'dBperkm':
        alpha_ph3 *= 434294.5

    return alpha_ph3
