import h2o_ddb
import h2o_bk
import matplotlib.pyplot as plt
import numpy as np

f = np.arange(1, 1000, 1)

P = 0.1
T = 300.0
X_h2o = 3.3E-04
X_he = 0.1
X_h2 = 1.0 - (X_h2o+X_he)
X = [X_h2, X_he, X_h2o]
P_dict = {'H2':0,'HE':1,'H2O':2}
otherpar = []
addb = h2o_ddb.alpha(f,T,P,X,P_dict,otherpar)
abk = h2o_bk.alpha(f,T,P,X,P_dict,otherpar)
#plt.loglog(f,addb)
plt.loglog(f,abk)
