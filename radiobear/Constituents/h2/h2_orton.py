import math
import os.path
import numpy as np
from scipy.interpolate import interp1d
import h2_jj


T0 = 273.0
atm2bar = 1.01325  # convert atm to bars
GHz = 29.9792458      # convert cm^-1 to GHz
tablename = 'orton_H2.tables'
h2tableList = {'eh2h2':0, 'nh2h2':1, 'eh2he':2, 'nh2he':3, 'eh2ch4':4, 'nh2ch4':5}
h2stateList = {'e':0, 'n':1}
Ttab = []
h2vab = []

def readInputFiles(freqs,path='./',verbose=False):
    ###===>overriding verbose!!!
    verbose = False
    global rpath
    rpath = path
    filename=os.path.join(rpath,tablename)
    ifp = open(filename,'r')

    # Temperature data read in as log save as linear
    data = ifp.readline().split()
    nTemp = int(data[0])
    Tmax = float(data[1])
    lTmx = math.log(Tmax)
    Tmin = float(data[2])
    lTmn = math.log(Tmin)
    dlT = (lTmx - lTmn)/(nTemp-1.0)
    ta = []
    ta.append(lTmn)
    for i in range(nTemp-1):
        ta.append(ta[i] + dlT)
    global Ttab
    for i in range(len(ta)):
        ta[i] = math.exp(ta[i])
    Ttab = np.array(ta)

    # frequency data
    nfreq=int(ifp.readline())
    data = ifp.readline().split()
    nperl = len(data)  #number of freq entries per line
    nfline = int(math.ceil(1.0*nfreq/nperl))  #number of lines
    fa = []
    for v in data:  #do first line and convert to GHz
        fa.append(float(v)*GHz)
    for i in range(nfline-1): #do the rest (and convert to GHz)
        data = ifp.readline().split()
        for v in data:
            fa.append(float(v)*GHz)
    ftab = np.array(fa)

    if verbose:
        print 'file = '+tablename
        print '\t%d frequencies per line and %d lines' % (nperl,nfline)
        print '\t%d frequencies between %.2f - %.2f GHz' % (nfreq,ftab[0],ftab[-1])
        print '\t%d temperatures between %.1f - %.1f K' % (nTemp,Tmin,Tmax)
        for i in range(10):
            print '\t\tftab[%d] = %f GHz' % (i,ftab[i])
        print '\t\t..... + %d others' % (nfreq-10)
        for i in range(len(Ttab)):
            print '\t\tTtab[%d] = %.4f' % (i,Ttab[i])
        print '\t\tEND'
    reallyVerbose=False

    # get temperature data for frequencies that are being used
    #     non-optimized approach rather than working monotonically through file
    #     but we only have to do it once per freq vector
    ms = (len(h2tableList), len(freqs), nTemp)
    global h2vab
    h2vab = np.zeros(ms)
    for ii in range(len(h2tableList)):
        for jj,f in enumerate(freqs):
            ifp.seek(0)
            ifreq = np.where(ftab>f)[0][0]
            if reallyVerbose:
                print 'f=%f  ftab=%f' % (f,ftab[ifreq])
            if ifreq==0: #if below lowest tabulated value extrapolate down
                ifreq=1
            nline2go =  2 + nfline + nfreq*ii + ifreq
            for kk in range(nline2go):
                v1 = ifp.readline()
            v1 = v1.split()
            v2 = ifp.readline().split() # we now have the bracketing values v1, v2
            v3 = ifp.readline().split() #  ... and now we can fit quadratic

            # fit quadratic, save as linear not log
            v = []
            X1   = ftab[ifreq-1]
            X21  = ftab[ifreq]   - ftab[ifreq-1]
            X32  = ftab[ifreq+1] - ftab[ifreq]
            X212 = ftab[ifreq]**2   - ftab[ifreq-1]**2
            X322 = ftab[ifreq+1]**2 - ftab[ifreq]**2
            for ll in range(len(v1)):
                Y1 = math.exp(float(v1[ll]))
                Y21 = (math.exp(float(v2[ll])) - Y1)
                Y32 = (math.exp(float(v3[ll])) - math.exp(float(v2[ll])))
                DQ = X212*X32 - X322*X21
                AQ = (X32*Y21 - X21*Y32)/DQ
                BQ = (X212*Y32 - X322*Y21)/DQ
                CQ = Y1 - AQ*X1**2 - BQ*X1
                v.append(AQ*f**2 + BQ*f + CQ)
            h2vab[ii,jj] = v
            if reallyVerbose:
                print 'nline2go = ',nline2go
                print 'v1:  ',v1
                print 'v2:  ',v2
                print 'v3:  ',v3
                print 'v:  ',v
    ifp.close()
    if reallyVerbose:
        print h2vab
    return 1

def alpha(freq,T,P,X,X_dict,otherPar,units='dBperkm',path='./',verbose=False):
    """Piece-wise quadratic frequency interpolation and extrapolation to orton_H2.tables.
       Interpolation in T is interptype ('cubic')
       Extrapolation in T is based on:
           down = nexp (=4) extrapolation in Temperature
           up = h2_jj matched at the nearest tabulated point.
               freq is a list of frequencies (must be a list)
               T and P are scalars in K and bars
               X is a list of the mixing ratios constituents at X_dict
               otherPar['h2state'] = 'e' or 'n' (h2states)
               Units either 'dBperkm' or 'invcm'
               Need to have otherPar['newset']=True whenever the frequencies change"""

    newset = otherPar['h2newset']
    h2state = otherPar['h2state']
    interpType = 'cubic'
    
    if len(Ttab)==0 or newset:
        readInputFiles(freq,path,verbose)

    h2state = h2state.lower()
    P_h2 = P*X[X_dict['H2']]
    P_he = P*X[X_dict['HE']]
    try:
        P_ch4 = P*X[X_dict['CH4']]
    except keyError:
        P_ch4 = 0.0

    alpha_h2 = []
    if T < Ttab[0]:        ### extrapolate down in T with quartic
        nearest = 0
        nexp = 4.0
        indi = [nearest,nearest+1]
        X1  = Ttab[indi[0]]
        X21n= Ttab[indi[1]]**nexp - Ttab[indi[0]]**nexp
        for ii in range(len(freq)):
            v = []
            for jj in indi:
                Tjj = Ttab[jj]
                #h2-h2
                xx = h2tableList['eh2h2'] + h2stateList[h2state]
                ah2near = h2vab[xx,ii,jj]
                #h2-he
                xx = h2tableList['eh2he'] + h2stateList[h2state]
                ahenear = h2vab[xx,ii,jj]
                #h2-ch4
                xx = h2tableList['eh2ch4'] + h2stateList[h2state]
                ach4near = h2vab[xx,ii,jj]
                #total
                a = (P_h2/atm2bar)*(ah2near*P_h2/atm2bar + ahenear*P_he/atm2bar + ach4near*P_ch4/atm2bar)*(T0/Tjj)**2
                v.append(a)
            Y1  = v[indi[0]]
            Y21 = v[indi[1]] - v[indi[0]]
            AQ = -1.0*Y21/X21n
            CQ = Y1 + AQ*(X1**nexp)
            a = CQ - AQ*(T**nexp)
            alpha_h2.append(a)   
    elif T > Ttab[-1]:     ### extrapolate up in T using Joiner
        nearest = -1
        Tnear = Ttab[nearest]
        jjnear = h2_jj.alpha(freq,Tnear,P,X,X_dict,otherPar)
        jj = h2_jj.alpha(freq,T,P,X,X_dict,otherPar)
        for ii in range(len(freq)):
            #h2-h2
            xx = h2tableList['eh2h2'] + h2stateList[h2state]
            ah2near = h2vab[xx,ii,nearest]
            #h2-he
            xx = h2tableList['eh2he'] + h2stateList[h2state]
            ahenear = h2vab[xx,ii,nearest]
            #h2-ch4
            xx = h2tableList['eh2ch4'] + h2stateList[h2state]
            ach4near = h2vab[xx,ii,nearest]
            #total
            anear = (P_h2/atm2bar)*(ah2near*P_h2/atm2bar + ahenear*P_he/atm2bar + ach4near*P_ch4/atm2bar)*(T0/Tnear)**2
            a = jj[ii]*(anear/jjnear[ii])
            alpha_h2.append(a)
    else:                     ### within range of Orton's tables
        for ii in range(len(freq)):
            #h2-h2
            xx = h2tableList['eh2h2'] + h2stateList[h2state]
            fv = interp1d(Ttab,h2vab[xx,ii],kind=interpType)
            ah2 = fv(T)
            #h2-he
            xx = h2tableList['eh2he'] + h2stateList[h2state]
            fv = interp1d(Ttab,h2vab[xx,ii],kind=interpType)
            ahe = fv(T)
            #h2-ch4
            xx = h2tableList['eh2ch4'] + h2stateList[h2state]
            fv = interp1d(Ttab,h2vab[xx,ii],kind=interpType)
            ach4 = fv(T)
            #total
            a = (P_h2/atm2bar)*(ah2*P_h2/atm2bar + ahe*P_he/atm2bar + ach4*P_ch4/atm2bar)*(T0/T)**2
            alpha_h2.append(a)
            
    if units=='dBperkm':
        for ii,a in enumerate(alpha_h2):
            alpha_h2[ii] = a*434294.5
    elif units != 'invcm':
        print 'Assuming units are in invcm...not '+units

        
    return alpha_h2
