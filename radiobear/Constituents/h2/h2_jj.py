# Some constants
coef = 3.9522E-14     # coefficient from
T0 = 273.0            # reference temperature in K
from radiobear.Constituents import parameters


def alpha(freq, T, P, X, P_dict, otherPar, **kwargs):
    par = parameters.setpar(kwargs)
    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_ch4 = P * X[P_dict['CH4']]
    th = T0 / T
    alpha_h2 = []
    for f in freq:
        f2 = f**2
        cf = coef * f2 * P_h2
        a = cf * (P_h2 * pow(th, 3.12) + 1.382 * P_he * pow(th, 2.24) + 9.322 * P_ch4 * pow(th, 3.34))
        if par.units == 'dBperkm':
            a *= 434294.5
        alpha_h2.append(a)
    return alpha_h2
