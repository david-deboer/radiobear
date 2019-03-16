import numpy as np
import math

def alpha(f,T,P,X,P_dict,units='dBperkm',path='./',verbose=False):
    """This is Karpowicz's model from:
       http://users.ece.gatech.edu/~psteffes/palpapers/karpowicz_data/water_model/karpowicz_h2o_model.m
       This function also requires the vvwlineshape function originally written by Jim Hoffman (with removal of df factor).
       Also added shift term as is done in Rosenkranz's work.
       Rosenkranz 1998 Radio Science vol 33, pp919-928
       NB (DDB): f may be a list, but T, P, etc should be scalars"""

    P_h2 = P*X[P_dict['H2']]
    P_he = P*X[P_dict['HE']]
    P_h2o= P*X[P_dict['H2O']]
    mbars_to_bars=0.001
    inv_km_to_dB=4.342945
    convert_to_km=1e-4
    To=300.0
    Theta=To/T
    NA=6.0221415e23             #molecules/mol
    M_amu=8.314472/0.46151805   #g/mol
    isotope_partition=0.997317
    M_amu_he=4.0026
    M_amu_h2=2.01594
        
    density_h2o=(M_amu*P_h2o)/(8.314472e-5*T)
    density_he =(M_amu_he*P_he)/(8.314310e-5*T)
    density_h2 =(M_amu_h2*P_h2)/(8.314472e-5*T)
    density_h2o=isotope_partition*(density_h2o/M_amu)*NA*(1.0/1e6)  #need in molecules/cc
        
    #Center Frequencies in GHZ
    f_o=[22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, 620.7008, 752.0332, 916.1712]
    #Line intensities
    I_o=[0.1314E-13, 0.2279E-11, 0.8058E-13, 0.2701E-11, 0.2444E-10, 0.2185E-11, 0.4637E-12, 0.2568E-10, 0.8392E-12, 0.3272E-11, 0.6676E-12, 0.1535E-08, \
         0.1711E-10, 0.1014E-08, 0.4238E-10]
    #Temperature coefficients
    E_o=[2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,3.597, 2.379, 2.852, 0.159, 2.391, 0.396, 1.441]
    #self broadening parameters converted to bars
    w_stmp=[0.01349, 0.01466, 0.01057, 0.01381, 0.01454, 0.009715, 0.00788, 0.01275, 0.00983, 0.01095, 0.01313, 0.01405, 0.011836, 0.01253, 0.01275]
    w_s = []
    for w in w_stmp:
        w_s.append(w/mbars_to_bars)
    x_s=[0.61, 0.85, 0.54, 0.74, 0.89, 0.62, 0.50, 0.67, 0.65, 0.64, 0.72, 1.0, 0.68, 0.84, 0.78]
    #foreign gas broadening parameters
    w_h2=[2.395, 2.4000, 2.395, 2.395, 2.390, 2.395, 2.395, 2.395, 2.395, 2.395, 2.395, 2.395, 2.395, 2.395, 2.395]
    w_he=[0.67, 0.71, 0.67,0.67, 0.63, 0.67, 0.67, 0.67, 0.67,0.67,0.67, 0.67, 0.67,0.67, 0.67]
    x_h2=[0.900, 0.950, 0.900, 0.900, 0.850, 0.900, 0.900, 0.900, 0.900, 0.900,0.900,0.900,0.900,0.900,0.900]
    x_he=[0.515, 0.490, 0.515, 0.490, 0.540, 0.515, 0.515,0.515, 0.515,0.515,0.515,0.515,0.515,0.515,0.515]
    SR=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    ####Calculate the line-broadening terms
    expo = []
    S = []
    for i,eo in enumerate(E_o):
        expo.append(eo*(1.0-Theta))
        S.append(I_o[i]*(Theta**2.5)*math.exp(expo[i]))
    #df aka gamma delta-nu change aka pressure broadening term
    df = []
    for i,w in enumerate(w_s):
        df.append( w*P_h2o*Theta**x_s[i] + w_h2[i]*P_h2*Theta**x_h2[i] + w_he[i]*P_he*Theta**x_he[i] )
    #calculate van-vleck line contribution...returns array of line contributions
    FSsum = vvwlinecontribution_modified(f,f_o,df,SR,S)
    line_contribution=inv_km_to_dB*convert_to_km*density_h2o*FSsum

    ####Calculate the continuum terms
    Cf_he=((1.0/mbars_to_bars)**2)*1.03562010226e-10            #  (dB/km)/((GHz x kPa)^2)->db/km((GHz bars)^2)  #eqn 10 Rosenkranz,1998
    Cf_h2=((1.0/mbars_to_bars)**2)*5.07722009423e-11
    Cs1=  4.36510480961e-07*pow(Theta,13.3619799812)          #equation 5 (correction applied from ^+6,^-6) Rosenkranz,1998
    Cs2=  2.10003048186e-26*pow(P_h2o/mbars_to_bars,6.76418487001)*pow(Theta,0.0435525417274)
    P_f=P_he+P_h2

    #Continuum Terms Foreign, and self 
    #Cf=((1/mbars_to_bars)^2) * 5.43e-10; %  (dB/km)/((GHz x kPa)^2)->db/km((GHz bars)^2)  %eqn 10 Rosenkranz,1998
    #Cs=((1/mbars_to_bars)^2)* 1.8e-8*Theta.^(4.5);    %equation 5 (correction applied from ^+6,^-6) Rosenkranz,1998

    #Foreign_Continuum=Cf.*P_f.*P_h2o.*f.^2*Theta^3; %foreign continuum term from Eqn 6, Rosenkranz
    #Self_Continuum=Cs*P_h2o^2.*f.^2*Theta^3; %self continuum term from Eqn 6, Rosenkranz

    farr = np.array(f)
    Foreign_Continuum_he=Cf_he*P_he*P_h2o*(farr**2)*pow(Theta,3.0)   #foreign continuum term from Eqn 6, Rosenkranz
    Foreign_Continuum_h2=Cf_h2*P_h2*P_h2o*(farr**2)*pow(Theta,3.0)
    Foreign_Continuum=Foreign_Continuum_he+Foreign_Continuum_h2
    Self_Continuum=Cs1*((P_h2o/mbars_to_bars)**2)*(farr**2.0)+Cs2*(farr**2.0)    #*power(Theta,3) #self continuum term from Eqn 6, Rosenkranz

    alpha_h2o=line_contribution+inv_km_to_dB*Foreign_Continuum+inv_km_to_dB*Self_Continuum
    #alpha_h2o = inv_km_to_dB*Self_Continuum
    #alpha_h2o = line_contribution
    if units!='dBperkm':
        alpha_h2o=alpha_h2o/434294.5

    return np.ndarray.tolist(np.ndarray.flatten(alpha_h2o))

def vvwlinecontribution_modified(f,fo,df,shift,S):
    """Taken from http://users.ece.gatech.edu/~psteffes/palpapers/karpowicz_data/water_model/vvwlineshape_modified.m
       This function determines the Van,Velck,Wiesskopf (VVW) lineshape over
       a range of freqencies which are user defined as is the frequency step
       the resulting vector is passed to this subroutine as f.
       fo is a vector containing the center frequencies (GHz) of the spectral 
       lines of the given molecule (ie PH3) as supplied by the "Submillimeter,
       Millimeter, and Microwave Spectral Catalog" aka Poynter-Pickett Line Catalog
       but in GHz, not in MHz as provided by the catalog.
       df is elsewhere called d_nu in the main program and is the line half-width at
       half-maximum or simply linewidth
    
       Bryan Karpowicz: Removed df term so fvvw is actually fvvw, and not
       df*fvvw. Also added cutoff at 750 GHz for consistency with Rosenkranz.
    
       NOTE: All frequencies are or have been converted to inverse cm!!!
       fvvw(f,fo,df)=(1/pi)*(f/fo)^2*(    df          +         df	   )
                                       --------------     --------------      -       base
                                       (fo-f-shift)^2 +df^2     (fo+f+shift)^2 +df^2
       F =                A      *(     B          +         C       )
       NB (DDB):  I am using np.matrix here to implement Bryans's matlab version.
       Return only summed arrays however;  matrices are only internal to this module"""

    n = np.size(f)
    m = np.size(fo)
    fmat = np.matrix(f)
    fomat = np.transpose(np.matrix(fo))
    dfmat = np.transpose(np.matrix(df))
    dfarr = np.array(df)
    shiftmat = np.transpose(np.matrix(shift))
    Smat = np.transpose(np.matrix(S))
    #f(1) f(2) f(3) f(4) ....f(n) n times where n is the number of frequency steps
    #f(1) f(2) f(3) f(4)				in the observation range                            
    # ...
    #m times where m is the number of spectral lines

    nones=np.matrix(np.ones(n))
    mones=np.transpose(np.matrix(np.ones(m)))
    f_matrix=mones*fmat
    fo_matrix=fomat*nones
    df_matrix=dfmat*nones
    shift_matrix=shiftmat*nones
    S_matrix = Smat*nones
    
    base=dfarr/(562500.0+dfarr**2)
    basemat = np.transpose(np.matrix(base))
    basemat = basemat*nones
    A=np.square(np.divide(f_matrix, fo_matrix))/math.pi
    B=np.divide(df_matrix, np.square(f_matrix-fo_matrix-shift_matrix) + np.square(df_matrix))
    C=np.divide(df_matrix, np.square(f_matrix+fo_matrix+shift_matrix) + np.square(df_matrix))
    cut_off=750.0

    #filter1=np.abs(f_matrix-fo_matrix-shift_matrix) <=cut_off
    #filter2=np.abs(f_matrix+fo_matrix+shift_matrix) <=cut_off

    F = np.multiply(A, B-basemat + C-basemat)
    FS = np.multiply(S_matrix,F)

    FSsum = np.array(np.transpose(np.sum(np.transpose(FS),1)))

    return FSsum

    #Note: the extra df factor comes in because of the expression for alpha has
    #an additional df factor which cancels with the alpha_max, but is added here
    #for computational reasons.  result is actually df*Fvvw

                                      
                                     
