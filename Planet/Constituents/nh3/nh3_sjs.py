from __future__ import absolute_import, division, print_function
import os.path
import numpy as np

# Some constants
coef = 7.244E+21     # coefficient from GEISA_PP.TEX eq. 14
T0 = 296.0           # reference temperature in K
hck = 1.438396       # hc/k  [K cm]
GHz = 29.9792458     # conversion from cm^-1 to GHz
fLower = 26.0        # beginning of transition between S and J (changed from 40 15/5/4)
fHigher = 34.0       # end of transition    "  (changed from 50 15/5/4)
EPS = 1E-12          # a small number

data = None


def readInputFiles(path, verbose=False):
    """This reads in the data files for nh3"""
    filename = os.path.join(path, 'nh3.npz')
    if verbose:
        print("Reading nh3 lines from {}".format(filename))
    global data
    data = np.load(filename)


def alpha(freq, T, P, X, P_dict, otherPar, units='dBperkm', path='./', verbose=False):
    Joiner = 0
    Spilker = 1
    Interp = 2

    # Read in data if needed
    if data is None:
        readInputFiles(path, verbose)

    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_nh3 = P * X[P_dict['NH3']]
    Pscale = 1.0 + P / 1.0E5  # ad hoc to make this match Berge-Gulkis for deep Jupiter atmosphere

    # Set Joiner
    GH2 = [1.690]
    GHe = [0.750]
    GNH3 = [0.6]
    ZH2 = [1.350]
    ZHe = [0.300]
    ZNH3 = [0.200]
    C = [1.0]
    D = [-0.45]
    # Set Spilker
    rexp = 8.79 * np.exp(-T / 83.0)
    GH2a = np.exp(9.024 - T / 20.3) - 0.9918 + P_h2
    try:
        GH2a = np.power(GH2a, rexp)
    except ValueError:
        GH2a = 0.0
    if GH2a < EPS:  # use Joiner regardless
        GH2.append(1.690)
        GHe.append(0.750)
        GNH3.append(0.60)
        ZH2.append(1.35)
        ZHe.append(0.30)
        ZNH3.append(0.20)
        C.append(1.00)
        D.append(-0.45)
    else:
        GH2a = 2.122 * np.exp(-T / 116.8) / GH2a
        GH2a = 2.34 * (1.0 - GH2a)
        GH2.append(GH2a)
        GHe.append(0.46 + T / 3000.0)
        GNH3.append(0.74)
        ZH2.append(5.7465 - 7.7644 * GH2a + 9.1931 * GH2a**2 - 5.6816 * GH2a**3 + 1.2307 * GH2a**4)
        ZHe.append(0.28 - T / 1750.0)
        ZNH3.append(0.50)
        C.append(-0.337 + T / 110.4 - T**2 / 70600.0)
        D.append(-0.45)
    # Set interp index
    GH2.append(0.0)
    GHe.append(0.0)
    GNH3.append(0.0)
    ZH2.append(0.0)
    ZHe.append(0.0)
    ZNH3.append(0.0)
    C.append(0.0)
    D.append(0.0)
    f0 = data['f0']
    I0 = data['I0']
    E = data['E']
    G0 = data['G0']
    n_dvl = 2.0 / 3.0
    n_int = 3.0 / 2.0
    ITG = I0 * np.exp(-((1.0 / T) - (1.0 / T0)) * E * hck)
    alpha_nh3 = []
    for f in freq:
        f2 = f**2
        if f <= fLower:
            use = Spilker
        elif f >= fHigher:
            use = Joiner
        else:
            use = Interp
            flfh = (fLower - fHigher) / (f - fLower)
            GH2[Interp] = GH2[Spilker] + (GH2[Spilker] - GH2[Joiner]) / flfh
            GHe[Interp] = GHe[Spilker] + (GHe[Spilker] - GHe[Joiner]) / flfh
            GNH3[Interp] = GNH3[Spilker] + (GNH3[Spilker] - GNH3[Joiner]) / flfh
            ZH2[Interp] = ZH2[Spilker] + (ZH2[Spilker] - ZH2[Joiner]) / flfh
            ZHe[Interp] = ZHe[Spilker] + (ZHe[Spilker] - ZHe[Joiner]) / flfh
            ZNH3[Interp] = ZNH3[Spilker] + (ZNH3[Spilker] - ZNH3[Joiner]) / flfh
            C[Interp] = C[Spilker] + (C[Spilker] - C[Joiner]) / flfh
            D[Interp] = D[Spilker] + (D[Spilker] - D[Joiner]) / flfh
        delta = D[use] * P_nh3

        gamma = np.power((T0 / T), n_dvl) * (GH2[use] * P_h2 + GHe[use] * P_he + G0 * GNH3[use] * P_nh3)
        g2 = gamma**2
        zeta = np.power((T0 / T), n_dvl) * (ZH2[use] * P_h2 + ZHe[use] * P_he + G0 * ZNH3[use] * P_nh3)
        z2 = zeta**2

        num = (gamma - zeta) * f2 + (gamma + zeta) * (np.power(f0 + delta, 2.0) + g2 - z2)
        den = np.power((f2 - np.power(f0 + delta, 2.0) - g2 + z2), 2.0) + 4.0 * f2 * g2
        shape = GHz * 2.0 * np.power(f / f0, 2.0) * num / (np.pi * den)
        alpha_nh3.append(np.sum(shape * ITG))

    alpha_nh3 = coef * (P_nh3 / T0) * pow((T0 / T), n_int + 2) * np.array(alpha_nh3) * Pscale
    if units == 'dBperkm':
        alpha_nh3 *= 434294.5

    return alpha_nh3
