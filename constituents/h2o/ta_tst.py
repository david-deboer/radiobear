import h2o_bk
import h2o_ddb
import matplotlib.pyplot as plt
import math
import numpy as np

jupiterTestP = [1.,10.,20.,100.,200.,400.,1000.,2000.,4000.,7500.,18000.]
jupiterTestT = [150.,300.,400.,700.,800.,1000.,1200.,1600.,2000.,2400.,3200.,]


otherPar = []


f = [4.0]#, 8.0,20.0]
fclr = ['b','r','g']


jupiterTestP = [1.,2.,4.,10.,20.,100.,200.,400.,1000.,2000.,4000.,7500.,18000.]
jupiterTestT = [150.,200.,250.,300.,400.,700.,800.,1000.,1200.,1600.,2000.,2400.,3200.,]
X_h2o = 0.005
X_he = 0.10
X_h2 = 1.0 - (X_h2o+X_he)
X_partial=[X_h2,X_he,X_h2o]
P_dict = {'H2':0,'HE':1,'H2O':2}

a_bk = []
a_ddb = []
for i in range(len(jupiterTestP)):
    T = jupiterTestT[i]
    P = jupiterTestP[i]
    a_bk.append(h2o_bk.alpha(f,T,P,X_partial,P_dict,otherPar))
    a_ddb.append(h2o_ddb.alpha(f,T,P,X_partial,P_dict,otherPar))
    print '%5.0f\t%4.0f\t%.3f\t%.3f' % (T,P,a_bk[i][0],a_ddb[i][0])
a_bk = np.array(a_bk)
a_ddb = np.array(a_ddb)
plt.figure('samples')
for i in range(len(f)):
    s = '%.0f GHz' % (f[i])
    print s
    plt.plot(jupiterTestP,a_bk[:,i],'bo')
    plt.plot(jupiterTestP,a_bk[:,i],'b',label=s)
    plt.plot(jupiterTestP,a_ddb[:,i],'ro')
    plt.plot(jupiterTestP,a_ddb[:,i],'r',label=s)
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('Pressure [bars]')
plt.ylabel(r'$\alpha$ [dB/km]')
s = r'%.1f GHz:  %.3f H$_2$, %.3f He, %.3f H$_2$O' % (f[0],X_h2,X_he,X_h2o)
plt.title(s)
print X_partial

if True:
    fp = open('juptest.dat','w')
    f = [1.4,4.0, 8.0,20.0]
    fclr = ['k','b','r','g']
    jupdat = np.loadtxt('jupiter.paulSolar')
    jupiterTestP = jupdat[:,2]
    jupiterTestT = jupdat[:,1]
    X_h2o = jupdat[:,7]
    X_he = jupdat[:,4]
    X_h2 = jupdat[:,3]
    a_bk = []
    a_ddb = []
    for i in range(len(jupiterTestP)):
        T = jupiterTestT[i]
        P = jupiterTestP[i]
        X_partial = [X_h2[i],X_he[i],X_h2o[i]]
        a_bk.append(h2o_bk.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_ddb.append(h2o_ddb.alpha(f,T,P,X_partial,P_dict,otherPar))
        fp.write('%f\t%f\t%f\t%f\t%f\t%s\t%s\n' % (P,T,X_partial[0],X_partial[1],X_partial[2],str(a_bk[i]),str(a_ddb[i])))
    a_bk = np.array(a_bk)
    a_ddb = np.array(a_ddb)
    plt.figure('ssampless')
    for i in range(len(f)):
        plt.plot(jupiterTestP,a_bk[:,i],fclr[i])
        plt.plot(jupiterTestP,a_ddb[:,i],fclr[i]+'--')
    plt.xscale('log')
    plt.xlabel('Pressure [bars]')
    plt.ylabel(r'$\alpha$ [dB/km]')
    #s = r'%.1f GHz:  %.3f H$_2$, %.3f He, %.3f NH$_3$' % (f[0],X_h2,X_he,X_nh3)
    #plt.title(s)
    print X_partial
    fp.close()
