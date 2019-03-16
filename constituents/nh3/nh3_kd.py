#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% THIS FUNCTION COMPUTES AMMONIA ABSORPTION UNDER JOVIAN ATMOSPHERIC
#% CONDITIONS BETWEEN 1-200 GHz
#%
#% NAME:
#%      NH3_Consistent_Model
#%
#% EXPLANATION:
#%   This function can be used to compute the opacity of pure ammonia and
#%   that of ammonia in a hydrogen/helium atmosphere between 0.1-200 GHz at
#%   pressures between 0.01-100 bars and temperatures between 200-500 K.
#%
#% CALLING SEQUENCE:
#%       alphanh3=NH3_Consistent_Model(f,T,P,H2mr,Hemr,NH3mr)
#%
#% INPUTS:
#%       f        - Array of frequencies (GHz)
#%       T        - Temperature (K)
#%       P        - Pressure (bars)
#%       H2mr     - Hydrogen mixing ratio (as mole fraction)
#%       Hemr     - Helium mixing ratio (as mole fraction)
#%       NH3mr    - Ammonia mixing ratio (as mole fraction)
#%
#% OUTPUTS:
#%      alphanh3  - Array of ammonia opacity (dB/km) at the input frequencies
#%
#% METHOD:
#%   A modified Ben Reuven (Ben Reuven, 1966) lineshape is used for computing
#%   ammonia opacity due to inversion lines, and a Gross lineshape (Gross, 1955)
#%   is used for computing ammonia opacity due to the rotational lines and
#%   the v2 roto-vibrational lines.
#%
#% History:
#%       written by Kiruthika Devaraj at Georgia Tech,  June, 2011
#%       converted to python by David DeBoer at UC Berkeley,  July 2012
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import math
import numpy as np
import os.path
import traceback

# Declare data arrays
fo = []
Io = []
Eo = []
gammaNH3o = []
H2HeBroad = []
fo_rot = []
Io_rot = []
Eo_rot = []
gNH3_rot = []
gH2_rot = []
gHe_rot = []
fo_v2 = []
Io_v2 = []
Eo_v2 = []

#%% Declaring Constants
GHztoinv_cm=1/29.9792458           #% for converting GHz to inverse cm
OpticaldepthstodB=434294.5	   #% convert from cm^-1 to dB/km
torrperatm=760.0                   #% convert from atm to torr
bartoatm=0.987                     #% convert from bat to atm
GHztoMHz=1000.0                    #% convert from GHz to MHz
hc=19.858252418E-24                #% planks (J.s) light (cm/s)
kB=1.38E-23                       #% boltzmann's in J/K or N.m/K
No=6.02297E23                      #% Avogadros Number [mole^-1]
R=8.31432E7                        #% Rydberg's [erg/mole-K]
To=300.0                           #% Ref temp for P/P Catalogue
dynesperbar=1.0E6                  #% dyne=bar/1e6;
coef=dynesperbar*No/R              #% See Appendix D: Using the Poyter-Pickett Catalogs

def viewsum(arr):
    """diagnostic viewer"""
    stack = traceback.extract_stack()
    filename, lineno, function_name, code = stack[-2]
    print '-------------',
    print code
    arrsum=0.0
    for a in arr:
        arrsum+=a
    print arrsum
    print np.shape(arr)
    print '+++++++++++++'

def readInputFiles(path,verbose=False):
    global fo, Io, Eo, gammaNH3o, H2HeBroad
    global fo_rot, Io_rot, Eo_rot, gNH3_rot, gH2_rot, gHe_rot
    global fo_v2, Io_v2, Eo_v2
    
    #%% Inversion lines: 
    #% fo is frequency in GHz, Io is line intensity in cm^-1/(molecule./cm^2), 
    #% Eo is lower state energy in cm^-1, gammaNH3o and H2HeBroad are self and 
    #% foreign gas broadening parameters.
    filename = os.path.join(path,'ammonia_inversion.dat')
    if verbose:
        print "Reading nh3 inversion lines:  "+filename
    fo,Io,Eo,gammaNH3o,H2HeBroad = np.loadtxt(filename,skiprows=1,unpack=True)
    nlin = len(fo)
    if verbose:
        print str(nlin)+' lines'
        
    #%% Rotational lines: 
    #% fo_rot is frequency in GHz, Io_rot is line intensity in 
    #% cm^-1/(molecule./cm^2), Eo_rot is lower state energy in cm^-1, gNH3_rot,
    #% gH2_rot, gHe_rot are broadening parameters for rotational lines.
    filename=os.path.join(path,'ammonia_rotational.dat')
    if verbose:
        print "Reading nh3 rotational lines:  "+filename
    fo_rot,Io_rot,Eo_rot,gNH3_rot,gH2_rot,gHe_rot = np.loadtxt(filename,skiprows=1,unpack=True)
    nlin_rot = len(fo_rot)
    if verbose:
        print str(nlin_rot)+' lines'
        
    #%% v2 roto-vibrational lines: 
    #% fo_v2 is frequency in GHz, Io_v2 is line intensity in
    #% cm^-1/(molecule./cm^2), Eo_v2 is lower state energy in cm^-1,
    filename = os.path.join(path,'ammonia_rotovibrational.dat')
    if verbose:
        print "Reading nh3 roto-vibrational lines:  "+filename
    fo_v2,Io_v2,Eo_v2 = np.loadtxt(filename,skiprows=1,unpack=True)
    nlin_v2 = len(fo_v2)
    if verbose:
        print str(nlin_v2)+' lines'
    return nlin,nlin_rot,nlin_v2

def alpha(freq,T,P,X,P_dict,otherPar,units='dBperkm',path='./',verbose=True):
    """function alphanh3=NH3_Consistent_Model(f,T,P,H2mr,Hemr,NH3mr)
    % The data files containing the frequency, line intensity and lower state 
    % energy for the ammonia transitions as given in the latest JPL spectral 
    % line catalog provided by Shanshan Yu and Brian Droiun  (personal 
    % communication, 2010), and the self and foreign gas broadening 
    % parameters for the various transitions as given by Devaraj 2011 are 
    % loaded in the following steps. """

    global fo, Io, Eo, gammaNH3o, H2HeBroad
    global fo_rot, Io_rot, Eo_rot, gNH3_rot, gH2_rot, gHe_rot
    global fo_v2, Io_v2, Eo_v2
    global GHztoinv_cm, OpticaldepthstodB, torrperatm, bartoatm
    global GHztoMHz, hc, kB, No, R, To, dynesperbar, coef

    if len(fo)==0:
        readInputFiles(path,verbose=verbose)

    P_h2 = P*X[P_dict['H2']]
    P_he = P*X[P_dict['HE']]
    P_nh3= P*X[P_dict['NH3']]

    #%% Computing partial pressures, temperature factor, and coefficient for ammonia
    #% Compute the mixing ratios of  of H2, He, and NH2
    H2mr = P_h2/P
    Hemr = P_he/P
    NH3mr = P_nh3/P
    #% Compute the temperature factor
    Tdiv=To/T
    #%  Coefficient for symmetric top molecule
    eta=3.0/2.0

    ####################INVERSION LINES
    #% Pressure Dependent Switch for the parameters of the inversion transitions
    #ddb 141218:  put in linear transition
    P_trans = 15.0
    dP_up = 5.0
    dP_down = 3.0
    if P>P_trans+dP_up:
        gnu_H2=1.6361;      gnu_He=0.4555;      gnu_NH3=0.7298;
        GAMMA_H2=0.8;       GAMMA_He=0.5;       GAMMA_NH3=1.0;
        zeta_H2=1.1313;     zeta_He=0.1;        zeta_NH3=0.5152;
        Z_H2=0.6234;        Z_He=0.5;           Z_NH3=2.0/3.0;
        d=0.2;
        Con=1.3746;
    elif P<=P_trans-dP_down:
        gnu_H2=1.7465;      gnu_He=0.9779;      gnu_NH3=0.7298;
        GAMMA_H2=0.8202;    GAMMA_He=1.0;       GAMMA_NH3=1.0;
        zeta_H2=1.2163;     zeta_He=0.0291;     zeta_NH3=0.5152;
        Z_H2=0.8873;        Z_He=0.8994;        Z_NH3=2.0/3.0;
        d=-0.0627;
        Con=0.9862;
    else: # in linear transition region
        gnu_H2   =1.7465 + (1.7465 - 1.6361)*(15.0-dP_down - P)/(dP_up + dP_down)
        gnu_He   =0.9779 + (0.9779 - 0.4555)*(15.0-dP_down - P)/(dP_up + dP_down)
        gnu_NH3  =0.7298
        GAMMA_H2 =0.8202 + (0.8202 - 0.8)*(15.0-dP_down - P)/(dP_up + dP_down)
        GAMMA_He =1.0    + (1.0    - 0.5)*(15.0-dP_down - P)/(dP_up + dP_down)
        GAMMA_NH3=1.0
        zeta_H2  =1.2163 + (1.2163 - 1.1313)*(15.0-dP_down - P)/(dP_up + dP_down)
        zeta_He  =0.0291 + (0.0291 -    0.1)*(15.0-dP_down - P)/(dP_up + dP_down)
        zeta_NH3 =0.5152
        Z_H2     =0.8873 + (0.8873 - 0.6234)*(15.0-dP_down - P)/(dP_up + dP_down)
        Z_He     =0.8994 + (0.8994 -    0.5)*(15.0-dP_down - P)/(dP_up + dP_down)
        Z_NH3    =2.0/3.0
        d        =-0.0627+ (-0.0627 -   0.2)*(15.0-dP_down - P)/(dP_up + dP_down)
        Con      =0.9862 + (0.9862 - 1.3746)*(15.0-dP_down - P)/(dP_up + dP_down)
        
    gammaNH3omat = np.matrix(gammaNH3o)
    #% Individual broadening parameters
    gH2=gnu_H2*P_h2
    gHe=gnu_He*P_he
    gNH3=gnu_NH3*P_nh3*gammaNH3omat
    #% Broadening parameter
    gamma=((gH2)*((Tdiv)**(GAMMA_H2))+(gHe)*((Tdiv)**(GAMMA_He))+gNH3*(295.0/T)**(GAMMA_NH3))
    #% Shift parameter
    delt=d*gamma
    #% Individual coupling parameters
    zH2=zeta_H2*P_h2
    zHe=zeta_He*P_he
    zNH3=zeta_NH3*P_nh3*gammaNH3omat
    #% Coupling parameter
    zeta=(zH2)*((Tdiv)**(Z_H2))+(zHe)*((Tdiv)**(Z_He))+zNH3*(295.0/T)**(Z_NH3)

    zetasize=np.size(fo)
    pst=delt      							#% answer in GHz
    #print delt
    #%Coupling element, pressure shift and dnu or gamma are in GHz, need to convert brlineshape to inverse cm which is done below
    n=np.size(freq)
    m=np.size(fo)
    fmat = np.matrix(freq)
    fomat = np.transpose(np.matrix(fo))
    #% f1 f2 f3 f4 ....fn  n times where n is the number of frequency steps
    #% f1 f2 f3 f4 ....fn				in the observation range
    #% ...
    #% f1 f2 f3 f4 ....fn
    #% m times where m is the number of spectral lines

    nones=np.matrix(np.ones(n))
    mones=np.transpose(np.matrix(np.ones(m)))
    f_matrix=mones*fmat
    fo_matrix=fomat*nones

    #% The 10^6 allows use of P(bar) for P(dynes/cm^2)
    expo = []
    ST = []
    alpha_noshape = []
    for i,eo in enumerate(Eo):
        expo.append(-(1.0/T-1.0/To)*eo*hc/kB)
        ST.append(Io[i]*math.exp(expo[i]))	                     #% S(T) =S(To)converted for temperature
        alpha_noshape.append(Con*coef*(P_nh3/To)*((To/T)**(eta+2.0))*ST[i])  #%0.9387
    #%Ben Reuven lineshape calculated by the brlineshape function gives the answer in GHz
    #%Here we change from GHz to inverse cm.
    
    dnu_matrix=np.transpose(gamma)*nones
    ce_matrix=np.transpose(zeta)*nones
    pst_matrix=np.transpose(pst)*nones
    Aa=(2.0/math.pi)*np.square(np.divide(f_matrix,fo_matrix))
    Bb=np.multiply((dnu_matrix-ce_matrix),np.square(f_matrix))
    Cc=dnu_matrix+ce_matrix
    Dd=np.square(fo_matrix+pst_matrix) + np.square(dnu_matrix)-np.square(ce_matrix)
    Ee=np.square(f_matrix)
    Jj=np.square(fo_matrix+pst_matrix)
    Gg=np.square(dnu_matrix)
    Hh=np.square(ce_matrix)
    Ii=4.0*np.multiply(np.square(f_matrix), np.square(dnu_matrix))
    Ff = np.divide(np.multiply(Aa,Bb+np.multiply(Cc,Dd)),np.square(Ee-Jj-Gg+Hh) + Ii)
    Fbr=(1.0/GHztoinv_cm)*Ff
    alpha_noshapemat = np.transpose(np.matrix(alpha_noshape))
    alpha_noshape_matrix=alpha_noshapemat*nones
    alpha_inversion=np.multiply(alpha_noshape_matrix,Fbr)
    
    ####################ROTATIONAL LINES
    ST_rot = []
    for i,eo in enumerate(Eo_rot):
        ST_rot.append(Io_rot[i]*(math.exp((1.0/To-1.0/T)*eo*hc/kB)))
    STmat_rot=np.transpose(np.matrix(ST_rot))
    #% Factor GAMMA:
    eta_H2=0.8730
    eta_He=2.0/3.0
    eta_NH3=1.0
    #% Factor nu:
    gnu_H2=0.2984
    gnu_He=0.75
    gnu_NH3=3.1789
    #% Factor Con
    Con=2.4268
    #% Individual broadening
    gH2=gnu_H2*P_h2*np.matrix(gH2_rot)
    gHe=gnu_He*P_he*np.matrix(gHe_rot)
    gNH3=gnu_NH3*P_nh3*np.matrix(gNH3_rot)
    #% Total broadening
    gamma_rot=((gH2)*((Tdiv)**(eta_H2))+(gHe)*((Tdiv)**(eta_He))+gNH3*(Tdiv)**(eta_NH3))
    n=np.size(freq)
    m=np.size(fo_rot)
    fmat = np.matrix(freq)
    fomat = np.transpose(np.matrix(fo_rot))
    nones=np.matrix(np.ones(n))
    mones=np.transpose(np.matrix(np.ones(m)))
    f_matrix=mones*fmat
    fo_matrix=fomat*nones
    dnu_matrix=np.transpose(gamma_rot)*nones;
    #% Gross Lineshape
    Aa=(4.0/math.pi)*np.multiply(np.square(f_matrix),dnu_matrix)
    Bb=np.square( np.square(fo_matrix)-np.square(f_matrix) )
    Cc=4.0*np.multiply(np.square(f_matrix),np.square(dnu_matrix))
    F_rot=np.divide(Aa,Bb+Cc);
    Fbr_rot=(1/GHztoinv_cm)*F_rot
    alpha_rot=np.multiply(Con*coef*(P_nh3/To)*((To/T)**(eta+2.0))*STmat_rot*nones,Fbr_rot)
    
    #%% Computing the opacity due to v2 roto-vibrational lines
    #% Computing the absorption contributed by the v2 rotovibrational lines
    ST_v2 = []
    for i,eo in enumerate(Eo_v2):
        ST_v2.append(Io_v2[i]*(math.exp((1./To-1./T)*eo*hc/kB)))
    STmat_v2=np.transpose(np.matrix(ST_v2))
    #% Broadening parameters for the v2 vibrational lines
    gH2_v2=np.transpose(np.matrix(np.linspace(1.4,1.4,len(fo_v2))))
    gHe_v2=np.transpose(np.matrix(np.linspace(0.68,0.68,len(fo_v2))))
    gNH3_v2=np.transpose(np.matrix(np.linspace(9.5,9.5,len(fo_v2))))
    #% Factor GAMMA
    eta_H2=0.73
    eta_He=0.5716
    eta_NH3=1.0
    #% Factor Con
    Con=1.1206
    #% Individual broadening parameters
    gH2=P_h2*gH2_v2
    gHe=P_he*gHe_v2
    gNH3=P_nh3*gNH3_v2
    #%Total broadening
    gamma_v2=((gH2)*((Tdiv)**(eta_H2))+(gHe)*((Tdiv)**(eta_He))+gNH3*(Tdiv)**(eta_NH3));
    n=np.size(freq)
    m=np.size(fo_v2)
    fmat = np.matrix(freq)
    fomat = np.transpose(np.matrix(fo_v2))
    nones=np.matrix(np.ones(n))
    mones=np.transpose(np.matrix(np.ones(m)))
    f_matrix=mones*fmat
    fo_matrix=fomat*nones
    dnu_matrix=gamma_v2*nones
    ##% Gross Lineshape
    Aa=(4.0/math.pi)*np.multiply(np.square(f_matrix),dnu_matrix)
    Bb=np.square( np.square(fo_matrix)-np.square(f_matrix) )
    Cc=4.0*np.multiply(np.square(f_matrix),np.square(dnu_matrix))
    F_v2=np.divide(Aa,Bb+Cc);
    Fbr_v2=(1/GHztoinv_cm)*F_v2
    alpha_v2=np.multiply(Con*coef*(P_nh3/To)*((To/T)**(eta+2.0))*STmat_v2*nones,Fbr_v2)
    
    #%% Computing the total opacity
    alpha_opdep=np.sum(np.transpose(alpha_inversion),1)+np.sum(np.transpose(alpha_rot),1)+np.sum(np.transpose(alpha_v2),1)
    alpha_opdep=np.array(np.transpose(alpha_opdep))
    if units == 'dBperkm':
        alpha_opdep*=OpticaldepthstodB
    alpha_nh3_temp=np.ndarray.tolist(np.ndarray.flatten(alpha_opdep))
    alpha_nh3 = []
    for aaa in alpha_nh3_temp:
        if aaa <= 0.0:
            aaa = 1E-8
        alpha_nh3.append(aaa)

    return alpha_nh3

    

