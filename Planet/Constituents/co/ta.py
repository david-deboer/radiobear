import co_ddb
import matplotlib.pyplot as plt
import numpy as np

#f = []
#for i in range(2000):
#    f.append(float(i+.01)/2.0)

f = np.linspace(115.2,115.4,1000)

P = 1.0
T = 300.0
X_co = 3.3E-5
X_he = 0.1
X_h2 = 1.0 - (X_co+X_he)
X = [X_h2, X_he, X_co]
P_dict = {'H2':0,'HE':1,'CO':2}
#a = co_ddb.alpha(f,T,P,X,P_dict,None)
#plt.semilogy(f,a,label='1')

#P=0.0009
#a = co_ddb.alpha(f,T,P,X,P_dict,None)
#plt.semilogy(f,a,label='9e-4')

P=0.0001
a = co_ddb.alpha(f,T,P,X,P_dict,None)
plt.semilogy(f,a,label='5e-2a')
##P=0.001
##avoigt = co_ddb.alpha(f,T,P,X,P_dict,'voigt')
##plt.semilogy(f,avoigt,label='5e-2voigt')
##P=0.001
##avvw = co_ddb.alpha(f,T,P,X,P_dict,'vvw')
##plt.semilogy(f,avvw,label='5e-2vvw')
##adiff = np.abs(np.array(avoigt) - np.array(avvw))
##plt.semilogy(f,adiff,label='diff')
##P=0.001
##a = co_ddb.alpha(f,T,P,X,P_dict,'diff')
##plt.semilogy(f,a,label='5e-2diff')
