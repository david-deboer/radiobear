from radiobear.constituents.nh3 import nh3_sjs
from radiobear.constituents.nh3 import nh3_hs
from numpy import exp


def alpha(freq, T, P, X, P_dict, other_dict, **kwargs):
    """
    This calls the appropriate formalism depending on the pressure.
    """
    PLower = 400.0
    PMid = 800.0
    PHigher = 2000.0
    wtype = 'linear'
    if P < PLower:
        alpha_nh3 = nh3_hs.alpha(freq, T, P, X, P_dict, other_dict, **kwargs)
    elif P > PHigher:
        alpha_nh3 = nh3_sjs.alpha(freq, T, P, X, P_dict, other_dict, **kwargs)
    else:
        a2 = nh3_sjs.alpha(freq, T, P, X, P_dict, other_dict, **kwargs)
        a1 = nh3_hs.alpha(freq, T, P, X, P_dict, other_dict, **kwargs)
        if wtype == 'exp':
            W = exp(-(P - PMid)**2 / 1000.0)
        else:
            W = (P - PLower) / (PHigher - PLower)
        alpha_nh3 = W * a2 + (1.0 - W) * a1
    return alpha_nh3
