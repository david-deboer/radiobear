from __future__ import print_function
from ta_jup_nh3_atm import *
import numpy as np
import matplotlib.pyplot as plt
import nh3_dbs
import nh3_bg
import nh3_dbs_sjs
import nh3_kd
import nh3_sjs
import nh3_sjsd

print('Available test functions:')
print('\n\tta.t0()\n\tta.t1()\n\tta.t2()\n\tta.t3()')

# 80% H2, 20% He, 1.5E-4 NH3, 1 bar 100K


def t0(fmin=1.0, fmax=300.0, fstep=1.0, x=[0.8, 0.2, 1.5e-4], P=1.0, T=100.0):
    P_dict = {'H2': 0, 'HE': 1, 'NH3': 2}
    otherPar = []
    f = np.arange(fmin, fmax + fstep, fstep)
    a = nh3_dbs_sjs.alpha(f, T, P, x, P_dict, otherPar)
    w = 30.0 / np.array(f)
    plt.loglog(w, a)


def t1(fmin=1.0, fmax=1000.0, fstep=1.0, Pmin=0.0001, Pmax=0.1, Pstep=0.001):
    f = np.arange(fmin, fmax + fstep, fstep)
    pv = np.arange(Pmin, Pmax, Pstep)
    for pressure_value in pv:
        P, T, X_partial = pressure_params(pressure_value)
        a = nh3_dbs_sjs.alpha(f, T, P, X_partial, P_dict, otherPar)
        plt.loglog(f, a)


def t2():
    f = [2.0, 8.0, 20.0, 40]
    ltyp = ['--', '-', '-', '--']

    a_kd = []
    a_bg = []
    a_sjs = []
    a_sjsd = []
    a_dbs = []
    a_dbs_sjs = []
    usenh3 = {'kd': [False, a_kd, nh3_kd.alpha, 'kd', 'k'],
              'bg': [True, a_bg, nh3_bg.alpha, 'b-g', 'r'],
              'sjs': [False, a_sjs, nh3_sjs.alpha, 'ts+jj/ps', 'b'],
              'sjsd': [False, a_sjsd, nh3_sjsd.alpha, 'kd+(ts+jj/ps)', 'g'],
              'dbs': [False, a_dbs, nh3_dbs.alpha, 'ab/ps', 'm'],
              'dbs_sjs': [True, a_dbs_sjs, nh3_dbs_sjs.alpha, 'ab/ps+(ts+jj/ps)', 'c']}

    V_Plt = []
    for i in range(len(P_Jup)):
        T = T_Jup[i]
        P = P_Jup[i]
        X_partial = X_Jup[i]
        V_Plt.append(P)
        for k in usenh3.keys():
            if usenh3[k][0]:
                usenh3[k][1].append(usenh3[k][2](f, T, P, X_partial, P_dict, otherPar))
    for k in usenh3.keys():
        if usenh3[k][0]:
            usenh3[k][1] = np.array(usenh3[k][1])
    plt.figure('samples')
    for k in usenh3.keys():
        if usenh3[k][0]:
            for i in range(len(f)):
                if i == 0:
                    s = usenh3[k][3]
                else:
                    s = None
                ll = usenh3[k][4] + ltyp[i]
                plt.plot(V_Plt, usenh3[k][1][:, i], ll, linewidth=2, label=s)
                Y = usenh3[k][1][10, i]
                s = '%.0f GHz' % (f[i])
                plt.text(V_Plt[10], Y, s)
                print(s)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.xlabel('Pressure [bars]')
    plt.ylabel(r'$\alpha$ [dB/km]')
    s = 'jupiter.paulvla10d'
    plt.title(s)


def t3():
    f = [1.4, 4.0, 8.0, 20.0, 40.0]
    fclr = ['c', 'b', 'r', 'g']
    a_kd = []
    a_bg = []
    a_sjs = []
    a_sjsd = []
    for i in range(len(P_Jup)):
        T = T_Jup[i]
        P = P_Jup[i]
        X_partial = X_Jup
        a_kd.append(nh3_kd.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_bg.append(nh3_bg.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_sjs.append(nh3_sjs.alpha(f,T,P,X_partial,P_dict,otherPar))
        a_sjsd.append(nh3_sjsd.alpha(f,T,P,X_partial,P_dict,otherPar))
        fp.write('%f\t%f\t%f\t%f\t%f\t%s\t%s\n' % (P,T,X_partial[0],X_partial[1],X_partial[2],str(a_kd[i]),str(a_bg[i])))
    a_kd = np.array(a_kd)
    a_bg = np.array(a_bg)
    a_sjs = np.array(a_sjs)
    a_sjsd = np.array(a_sjsd)
    plt.figure('ssampless')
    for i in range(len(f)):
        plt.plot(P_Jup,a_kd[:,i],fclr[i])
        plt.plot(P_Jup,a_bg[:,i],fclr[i]+'--')
        plt.plot(P_Jup,a_sjs[:,i],fclr[i]+':')
        plt.plot(P_Jup,a_sjs[:,i],'k')
    plt.xscale('log')
    plt.xlabel('Pressure [bars]')
    plt.ylabel(r'$\alpha$ [dB/km]')
    fp.close()
