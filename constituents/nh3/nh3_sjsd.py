import nh3_sjs, nh3_kd
import numpy as np

def alpha(freq,T,P,X,P_dict,otherPar,units='dBperkm',path='./',verbose=False):
    PLower = 10.0
    PMid = 35.0
    PHigher = 100.0
    wtype = 'linear'
    if P < PLower or P > PHigher:
        alpha_nh3 = nh3_sjs.alpha(freq,T,P,X,P_dict,otherPar,units,path,verbose)
    else:
        a1 = np.array(nh3_sjs.alpha(freq,T,P,X,P_dict,otherPar,units,path,verbose))
        a2 = np.array(nh3_kd.alpha(freq,T,P,X,P_dict,otherPar,units,path,verbose))
        if wtype=='exp':
            W = np.exp(-(P-PMid)**2/1000.0)
        else:
            if P < PMid:
                W = (P-PLower) / (PMid-PLower)
            else:
                W = 1 - (P-PMid) / (PHigher-PMid)
        alpha_nh3 = W*a2 + (1.0-W)*a1
    return alpha_nh3
