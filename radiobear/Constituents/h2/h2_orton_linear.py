import math
import os.path
import numpy as np
#import h2_jj

T0 = 273.0
atm2bar = 1.01325  # convert atm to bars
GHz = 29.9792458      # convert cm^-1 to GHz
tablename = 'orton_H2.tables'
h2etc = {'eh2h2':0, 'nh2h2':1, 'eh2he':2, 'nh2he':3, 'eh2ch4':4, 'nh2ch4':5}
h2type = {'e':0, 'n':1}
Ttab = []
h2vab = []

def readInputFiles(freqs,path='./',verbose=False):
    global rpath
    rpath = path
    filename=os.path.join(rpath,tablename)
    ifp = open(filename,'r')

    # Temperature data
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
            print '\t\tftab[%d]=%f GHz' % (i,ftab[i])
        for t in Ttab:
            print '\t\t%.4f ==> %.4f' % (t, math.exp(t))
    reallyVerbose=False

    # get temperature data for frequencies that are being used
    #     non-optimized approach rather than working monotonically through file
    #     but we only have to do it once per freq vector
    ms = (len(h2etc), len(freqs), nTemp)
    global h2vab
    h2vab = np.zeros(ms)
    for ii in range(len(h2etc)):
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

            v = []
            dfr = (f - ftab[ifreq-1])/(ftab[ifreq]-ftab[ifreq-1])
            for ll in range(len(v1)):
                dy = (float(v2[ll]) - float(v1[ll]))
                v.append(float(v1[ll]) + dy*dfr)
            h2vab[ii,jj] = v
            if reallyVerbose:
                print 'nline2go = ',nline2go
                print 'dfr: ',dfr
                print 'v1: ',v1
                print 'v2: ',v2
                print 'v: ',v
    ifp.close()
    if reallyVerbose:
        print h2vab
    return 1

def alpha(freq,T,P,X,X_dict,IQ,units='dBperkm',path='./',newset=True,verbose=False):
    """freq is a list of frequencies (must be a list)
       T and P are scalars in K and bars
       X is a list of the mixing ratios constituents at X_dict
       IQ = 'e' or 'n'
       Units either 'dBperkm' or 'invcm'
       Need to have reset=True whenever the frequencies change"""
    
    if len(Ttab)==0 or newset:
        readInputFiles(freq,path,verbose)

    IQ = IQ.lower()
    P_h2 = P*X[X_dict['H2']]
    P_he = P*X[X_dict['HE']]
    P_ch4 = P*X[X_dict['CH4']]

    alpha_h2 = []
    logT = math.log(T)
    for ii in range(len(freq)):
        #h2-h2
        xx = h2etc['eh2h2'] + h2type[IQ]
        v = np.interp(logT,Ttab,h2vab[xx,ii])
        ah2 = math.exp(v)
        #h2-he
        xx = h2etc['eh2he'] + h2type[IQ]
        v = np.interp(logT,Ttab,h2vab[xx,ii])
        ahe = math.exp(v)
        #h2-ch4
        xx = h2etc['eh2ch4'] + h2type[IQ]
        v = np.interp(logT,Ttab,h2vab[xx,ii])
        ach4 = math.exp(v)
        #total
        a = (P_h2/atm2bar)*(ah2*P_h2/atm2bar + ahe*P_he/atm2bar + ach4*P_ch4/atm2bar)*(T0/T)**2
        if units=='dBperkm':
            a*=434294.5
        elif units != 'invcm':
            print 'Assuming units are in invcm...not '+units
        alpha_h2.append(a)
        
    return alpha_h2
