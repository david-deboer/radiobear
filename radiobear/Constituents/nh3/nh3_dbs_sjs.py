from radiobear.Constituents.nh3 import nh3_sjs
from radiobear.Constituents.nh3 import nh3_dbs
from radiobear.Constituents import parameters


def alpha(freq, T, P, X, P_dict, otherPar, **kwargs):
    """
    This calls the appropriate formalism depending on the pressure.
    """
    PLower = 400.0
    PMid = 800.0
    PHigher = 2000.0
    wtype = 'linear'
    if P < PLower:
        alpha_nh3 = nh3_dbs.alpha(freq, T, P, X, P_dict, otherPar, **kwargs)
    elif P > PHigher:
        alpha_nh3 = nh3_sjs.alpha(freq, T, P, X, P_dict, otherPar, **kwargs)
    else:
        a2 = nh3_sjs.alpha(freq, T, P, X, P_dict, otherPar, **kwargs)
        a1 = nh3_dbs.alpha(freq, T, P, X, P_dict, otherPar, **kwargs)
        if wtype == 'exp':
            W = np.exp(-(P - PMid)**2 / 1000.0)
        else:
            W = (P - PLower) / (PHigher - PLower)
        alpha_nh3 = W * a2 + (1.0 - W) * a1
    return alpha_nh3
