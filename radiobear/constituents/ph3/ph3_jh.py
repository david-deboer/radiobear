from __future__ import absolute_import, division, print_function
import os.path
import numpy as np
from radiobear.constituents import parameters

# Some constants
coef = 7.244E+21        # coefficient from GEISA_PP.TEX eq. 14
T0 = 300.0              # reference temperature in K
hck = 1.438396          # hc/k  [K cm]
GHz = 29.9792458        # conversion from cm^-1 to GHz
PI = np.pi

# Set data arrays and values
data = None
data_wgt = {}


def readInputFiles(par):
    # Read in lines
    filename = os.path.join(par.path, 'ph3jh.npz')
    if par.verbose:
        print("Reading ph3 lines from  {}".format(filename))
    global data
    data = {}
    data_in = np.load(filename)
    for x in data_in.files:
        data[x] = data_in[x]
    if par.truncate_strength is not None:
        if par.verbose:
            print("Truncating lines less than {}".format(par.truncate_strength))
        used_I0 = np.where(data['I0'] > par.truncate_strength)
        data['f0'] = data['f0'][used_I0]
        data['I0'] = data['I0'][used_I0]
        data['E'] = data['E'][used_I0]
    test_length = len(data['f0'])
    if par.truncate_freq is not None:
        if par.verbose:
            print("Truncating lines greater than {}".format(par.truncate_freq))
        used_f = np.where(data['f0'] < par.truncate_freq)
        data['f0'] = data['f0'][used_f]
        data['I0'] = data['I0'][used_f]
        data['E'] = data['E'][used_f]
    # Read in weights
    filename = os.path.join(par.path, 'PH3WGT.npz')
    if par.verbose:
        print("Reading ph3 wgts from {}".format(filename))
    global data_wgt
    data_in = np.load(filename)
    for x in data_in.files:
        data_wgt[x] = data_in[x]
    if par.truncate_strength:
        data_wgt['WgtI0'] = data_wgt['WgtI0'][used_I0]
        data_wgt['WgtFGB'] = data_wgt['WgtFGB'][used_I0]
        data_wgt['WgtSB'] = data_wgt['WgtSB'][used_I0]
    if len(data_wgt['WgtI0']) != test_length:
        raise ValueError("PH3 wgt vector wrong length.")
    if par.truncate_freq:
        data_wgt['WgtI0'] = data_wgt['WgtI0'][used_f]
        data_wgt['WgtFGB'] = data_wgt['WgtFGB'][used_f]
        data_wgt['WgtSB'] = data_wgt['WgtSB'][used_f]
    del data_in


def alpha(freq, T, P, X, P_dict, other_dict, **kwargs):
    """Computes the absorption due to ph3"""

    par = parameters.setpar(kwargs)
    # Read in data if needed
    if data is None:
        readInputFiles(par)

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
    gamma = np.power((T0 / T), n_dvl) * (GH2 * P_h2 + GHe * P_he) * WgtFGB + np.power((T0 / T), 1.0) * GPH3 * P_ph3 * WgtSB  # noqa
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
    if par.units == 'dBperkm':
        alpha_ph3 *= 434294.5

    del num, den, shape
    return alpha_ph3
