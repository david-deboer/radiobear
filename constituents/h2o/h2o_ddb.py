import math
import os.path

# Some constants
T0 = 300.0           # reference temperature in K
AMU_H2O=        18.015
R=              8.314462E7

#Set data arrays
f0 = []
Ei = []
A = []
GH2 = []
GHe = []
GH2O = []
x_H2 = []
x_He = []
x_H2O = []

def readInputFiles(path,verbose=False):
    """If needed this reads in the data files for h2o"""
    useLinesUpTo = 10   # index number
    global nlin
    nlin = 0
    if verbose:
        print "Reading h2o lines"
    filename = os.path.join(path,'h2od.lin')
    ifp = open(filename,'r')
    for line in ifp:
        if nlin >= useLinesUpTo:
            break
        nlin+=1
        data = line.split()
        if len(data) == 9:
            f0.append(float(data[0]))
            Ei.append(float(data[1]))
            A.append(float(data[2]))
            GH2.append(float(data[3]))
            GHe.append(float(data[4]))
            GH2O.append(float(data[5]))
            x_H2.append(float(data[6]))
            x_He.append(float(data[7]))
            x_H2O.append(float(data[8]))
        else:
            break
    ifp.close()
    if verbose:
        print '   '+str(nlin)+' lines'
    return nlin

def alpha(freq,T,P,X,P_dict,otherPar,units='dBperkm',path='./',verbose=False):

    # Read in data if needed
    if len(f0)==0:
        readInputFiles(path,verbose)

    P_h2 = P*X[P_dict['H2']]
    P_he = P*X[P_dict['HE']]
    P_h2o= P*X[P_dict['H2O']]
    n_int = 3.0/2.0
    rho = 1.0E12*AMU_H2O*P_h2o/(R*T)
    Pa = 0.81*P_h2 + 0.35*P_he

    alpha_h2o = []
    for f in freq:
        f2 = f**2
        alpha = 0.0
        for i in range(nlin):
            gamma = pow((T0/T),x_H2[i])*GH2[i]*P_h2
            gamma+= pow((T0/T),x_He[i])*GHe[i]*P_he
            gamma+= pow((T0/T),x_H2O[i])*GH2O[i]*P_h2o
            g2 = gamma**2
            ITG = A[i]*math.exp(-Ei[i]/T)
            shape = gamma/( (f0[i]**2 - f2)**2 + 4.0*f2*g2)
            alpha += shape*ITG
        GR1971 = 1.08E-11*rho*pow((T0/T),2.1)*Pa*f2
        a = 2.0*f2*rho*pow((T0/T),n_int)*alpha/434294.5 + GR1971/434294.5
        if units=='dBperkm':
            a*=434294.5
        alpha_h2o.append(a)
        
    return alpha_h2o

