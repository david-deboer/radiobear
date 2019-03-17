# Some constants
coef = 3.9522E-14     # coefficient from
T0 = 273.0            # reference temperature in K


def alpha(freq, T, P, X, P_dict, otherPar, units='dBperkm', path='./', verbose=False):

    P_h2 = P * X[P_dict['H2']]
    P_he = P * X[P_dict['HE']]
    P_ch4 = P * X[P_dict['CH4']]
    th = T0 / T
    alpha_h2 = []
    for f in freq:
        f2 = f**2
        if otherPar['h2state'] == 'e':
            pre = (T / 55.0)**2.7
            if pre > 1.0:
                pre = 1.0
            pre *= (T / 120.0)**0.55
            if pre > 1.0:
                pre = 1.0
        elif otherPar['h2state'] == 'n':
            pre = (T / 40.0)**2.5
            if pre > 1.0:
                pre = 1.0
        else:
            print('INVALID H2STATE')
            return 0.0
        cf = coef * f2 * P_h2 * pre
        a = cf * (P_h2 * pow(th, 3.12) + 1.382 * P_he * pow(th, 2.24) + 9.322 * P_ch4 * pow(th, 3.34))
        if units == 'dBperkm':
            a *= 434294.5
        alpha_h2.append(a)

    return alpha_h2
