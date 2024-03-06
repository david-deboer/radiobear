# -*- mode: python; coding: utf-8 -*-
"""Modify atmosphere."""
# Copyright 2018 David DeBoer
# Licensed under the 2-clause BSD license.
import sys
import os
import numpy as np 
import matplotlib.pyplot as plt
from . import chemistry 


def tweakAtm(self):
    """Tweak the atmosphere data."""
    # Import tweakmodule
    __import__(self.config.tweakmodule)
    tweakModule = sys.modules[self.config.tweakmodule]

    # Run module then log
    self.tweakComment, self.gas, self.cloud = tweakModule.modify(self.gas, self.cloud,
                                                                 self.config.C, self.config.Cl)
    if self.verbose:
        print('---tweakComment')
        print(self.tweakComment)
        print('---')
    self.log.add(self.tweakComment, False)
    tf = os.path.join(self.config.path, self.config.tweakmodule + '.py')
    tp = open(tf, 'r')
    dt = tp.read()
    self.log.add('======================' + tf + '=====================', False)
    self.log.add(dt, False)
    self.log.add('====================================================================', False)
    tp.close()


def scaleAtm(self, scale_info='Scratch/scale.dat'):
    """Built-in tweak module."""
    if isinstance(scale_info, str):
        import alpha
        col, scale_info = alpha.read_scalefile(scale_info)
    else:
        col = scale_info.keys()

    if len(scale_info[col[0]]) != self.nAtm:
        print("Warning - scale file doesn't match atmosphere.  Not applying.")
        return None

    for i in range(self.nAtm):
        for gas in col:
            if gas.lower() != 'p':
                self.gas[self.config.C[gas.upper()]][i] *= scale_info[gas][i]


# From Asplund Solar abundance ratios to convert between solar enhancement and ppm mix
N2H = 1.48e-4 # N/H2 ratio 
S2H = 2.89e-5 # S/H2 ratio 
O2H = 1.07e-3 # O/H2 ratio 

def model_stochastic(self, mix=[0.3,-0.3,0.3,-0.3,0.3,-0.3,0.3,-1], P_n=[100,40,20,15,10,6,4,2,0.6], rh_h2o=1.0, rh_nh3=0.42, deep_nh3 = 2.30, deep_h2s = 2.30, deep_h2o=2.30, nh3_h2o = 0.03, mixing_mode='uniform', adiabat='dry', plotting=False ): 
    ''' Vertical atmosphere structure based on the thermo-chemical equlibrium calculations

    Parameters
    ----------
    planet : class 
        [-] Instance of a planet created in radiobear 
    Keyword Arguments
    ----------
    
    mix : nx1 float  
        [solar] Mixing ratio slope per log pressure 
    P_n : nx1 float   
        [bar] Pressures nodes for mixing  ratio, can extend higher than ammonia cloud condensation pressure 
    rh_nh3 : float 
        [-] Relative humidity for the condenstion of ammonia 
    rh_h2o : float 
        [-] Relative humidity for the condenstion of water 
    deep_nh3 : float 
        [solar] Deep ammonia abundance (corresponding to deepest pressure level)
    deep_h2s : float 
        [solar] Deep hydrogen sulfide abundance 
    deep_h2o: float 
        [solar] Deep water abundance 
    nh3_h2o : float
        [-]  Ammonia solubility in water 

    mixing_mode : str
        'coupled' : all traces gases have the same mixing gradients 
        'constant' : only nh3 is modified, while h2s and h2o are constant until condensation pressure 

    Returns
    -------
    Locally modifies the vertical abundance profiles 

    
    Warnings
    -------
    
    
    Example
    -------
    
    import radiobear as rb  
    import numpy as np 

    j = rb.planet.Planet('jupiter', plot_atm=False, plot_bright=False,verbose=False) 

    j.atmos[0].cloud[j.atmos[0].config.Cl['SOLN']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['SOLN']])
    j.atmos[0].cloud[j.atmos[0].config.Cl['H2O']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['H2O']])
    j.atmos[0].cloud[j.atmos[0].config.Cl['NH4SH']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['NH4SH']])
    j.atmos[0].cloud[j.atmos[0].config.Cl['NH3']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['NH3']])

    rb.atm_modify.model_stochastic(j, plotting=True,mixing_mode='constant')




    References
    ------------
    Moeckel et al., 2022 - Ammonia Abundance Derived from VLA and Juno observations 

    Todo
    -----
    Implement feedback loop between condensation and temperature structure  

    Notes
    -------
    12/16/19, CM, initial comit 
    4/2/21 CM, updated H2S reference value to be consistent with Solar abundances
    7/12/22 CM, Rewrote the function to fit within radiobear 
    '''
    # Check that input is consistent 
    assert len(mix)==len(P_n)-1 
    # Check that pressure is ascending 
    assert np.all(np.diff(P_n) < 0 )
    # Check that pressure nodes are sorted  

    P = self.atmos[0].gas[self.atmos[0].config.C['P']]

    T = self.atmos[0].gas[self.atmos[0].config.C['T']]

    ''' 
    To implement later 
    # # Import the difference between dry and wet adiabat 
    # if adiabat.lower() == 'dry':
    #     T = TPprofile(P,adiabat='dry', updated=True)
    # elif adiabat.lower() == 'wet':
    #     T = TPprofile(P,adiabat='wet', updated=True)
    # elif adiabat>0 and adiabat < 1.5: 

    #     Td  = TPprofile(adiabat='dry') 
    #     Tw  = TPprofile(adiabat='wet') 
    #     T = Tw + (Td-Tw)*adiabat 
    # else:
    #     try: 
    #         T = j.atmos[0].gas[j.atmos[0].config.C['T']]
    ''' 

    # 

    # Ammonia 
    # Compute the value under the cloud to keep constant between highest pressure node and condensation point 
    subcloud_nh3 = (deep_nh3 - np.dot(mix,-np.diff(np.log(P_n))))
    nh3_prof = np.ones_like(P)*subcloud_nh3*N2H

    # Hydrogn sulfide 
    # Scale thermoequilbrium profile to deep_h2s. H2S must be taking out since it's never been observed 
    h2s_prof = np.copy(self.atmos[0].gas[self.atmos[0].config.C['H2S']])*(deep_h2s*S2H)/self.atmos[0].gas[self.atmos[0].config.C['H2S']][-1]
    
    # Water 
    # If mixing_mode=constant, h2o profile is set to deep water 
    # If mixing_mode=uniform, h2o profile follows the mixing regime of ammonia 
    h2o_prof = np.ones_like(P)*deep_h2o*O2H


    # Set the ammonia abundance in the stratosphere to zero 
    nh3_prof[P<0.2] = 1e-13
    h2s_prof[P<0.2] = 1e-13
    h2o_prof[P<0.2] = 1e-13

    # Copmute the deep atmosphere until you hit the last pressure node  
    n_nodes = len(P_n)

    idx_n = np.zeros(n_nodes)
    idx_n[0] =  int(np.where(P>P_n[0])[0][0]) # first index for the pressure node 
    # Set the deep abundances 
    nh3_prof[int(idx_n[0]):] = deep_nh3*N2H
    for i in range(1,n_nodes): 
        idx_n[i] = int(np.where(P>P_n[i])[0][0])
        nh3_prof[int(idx_n[i]):int(idx_n[i-1])] = nh3_prof[int(idx_n[i-1])] - (np.log(P_n[i-1]) - np.log(P[int(idx_n[i]):int(idx_n[i-1])]))*mix[i-1]*N2H 
        if mixing_mode.lower().strip() == 'coupled': 
            #h2s_prof[int(idx_n[i]):int(idx_n[i-1])] = h2s_prof[int(idx_n[i-1])] - (np.log(P_n[i-1]) - np.log(P[int(idx_n[i]):int(idx_n[i-1])]))*mix[i-1]*S2H  
            h2o_prof[int(idx_n[i]):int(idx_n[i-1])] = h2o_prof[int(idx_n[i-1])] - (np.log(P_n[i-1]) - np.log(P[int(idx_n[i]):int(idx_n[i-1])]))*mix[i-1]*O2H 


    # Loop between 0.2 and 10 bar for condensation (Note this only works for Jupiter)
    # i_s =  int(np.where(P>0.2)[0][0]) # first index for the pressure node 
    # i_f =  int(np.where(P>10)[0][0]) # last index for the pressure node 

    for i in range(len(P)):
        # Calculate the saturation pressure for the relevant gases at the given pressure 
        # --------------------------------------------------------------
        # NH3 
        psat_nh3 = chemistry.ConstituentProperties('NH3').Psat(T[i])
        pp_nh3= P[i] * nh3_prof[i]

        # H2S 
        psat_h2s = chemistry.ConstituentProperties('H2S').Psat(T[i])
        pp_h2s= P[i] * h2s_prof[i]

        # H2O profile 
        psat_h2o = chemistry.ConstituentProperties('H2O').Psat(T[i])
        pp_h2o = P[i] * h2o_prof[i]


        # Make the NH3 clouds 
        # --------------------------------------------------------------
        cpp_nh3 = nh3_prof[i]*P[i]

        # The relative humidity determines if over or undersaturate, and then changes the saturation pressure accordingly 
        if cpp_nh3 > psat_nh3*rh_nh3:
            nh3_prof[i] = rh_nh3 * psat_nh3/P[i]

        # Make the H2S clouds 
        # --------------------------------------------------------------
        cpp_h2s = h2s_prof[i]*P[i]

        # The relative humidity determines if over or undersaturate, and then changes the saturation pressure accordingly 
        if cpp_h2s > psat_h2s:
            h2s_prof[i] = psat_h2s/P[i]

        # Make the H2O clouds 
        # --------------------------------------------------------------
        cpp_h2o = h2o_prof[i]*P[i]

        # The relative humidity determines if over or undersaturate, and then changes the saturation pressure accordingly 
        if cpp_h2o > psat_h2o*rh_h2o:
            h2o_prof[i] = rh_h2o*psat_h2o/P[i]

    if plotting: 
        fig,ax = plt.subplots(figsize=(9,9)) 
        ax.plot(self.atmos[0].gas[self.atmos[0].config.C['NH3']]/N2H,P,color='gray', label=r'NH$_3$')
        ax.plot(self.atmos[0].gas[self.atmos[0].config.C['H2S']]/S2H,P,color='yellow', label=r'H$_2$S')
        ax.plot(self.atmos[0].gas[self.atmos[0].config.C['H2O']]/O2H,P,color='navy', label=r'H$_2$O')    


        ax.plot(nh3_prof/N2H,P,linestyle='--',color='gray', label=r'Modified NH$_3$')
        ax.plot(h2s_prof/S2H,P,linestyle='--',color='yellow', label=r'Modified H$_2$S')
        ax.plot(h2o_prof/O2H,P,linestyle='--',color='navy', label=r'Modified H$_2$O')

        ax.set_yscale('log')
        ax.invert_yaxis()
        ax.set_ylim([P_n[0],0.1])

        ax.set_ylabel('P [bar]')
        ax.set_xlabel('Abundance [solar]')
        plt.legend()
        plt.show()

    # Updates the ammonia abundance 
    self.atmos[0].gas[self.atmos[0].config.C['NH3']] = np.copy(nh3_prof)
    self.atmos[0].gas[self.atmos[0].config.C['H2S']] = np.copy(h2s_prof)
    self.atmos[0].gas[self.atmos[0].config.C['T']]   = np.copy(T)
    self.atmos[0].gas[self.atmos[0].config.C['H2O']] = np.copy(h2o_prof)

    return None


#def nh3_profile(j, deep_nh3=340.56e-6, deep_h2s=None, deep_h2o=None, rh_nh3 = 0.42, rh_nh4sh = 0.08, nh3_h2o = 0.03, dnh3dlP= 50e-6 , Pt_vm=15.5, Pb_vm=50,  adiabat='dry',  plotting=False, scalefactor=False, figsize=(5,5)):

def model_tceq(self,deep_nh3 = 2.30, deep_h2s = 2.30, deep_h2o = 2.30, rh_nh3 = 0.42, rh_nh4sh = 0.08, nh3_h2o = 0.03, rh_h2o = 1.0, dnh3dlP = 0.3, Pt_vm = 15.5, Pb_vm = 30, P_stra=0.2, plotting = False, adiabat = 'dry'):
    ''' 
    # Atmospheric structure based on modified thermo-chemical equilibrium models 


    Parameters
    ----------
    j : x class 
        [-] Radio bear instance of a planet 

    Keyword Arguments
    ----------

    list all here 

    Returns
    -------
    Modifys the ammonia, h2s profile in the planet instance j 

    
    Warnings
    -------
    
    
    Example
    -------
    import numpy as np 
    import radiobear as rb  

    j = rb.planet.Planet('jupiter', plot_atm=False, plot_bright=False,verbose=False) 

    j.atmos[0].cloud[j.atmos[0].config.Cl['SOLN']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['SOLN']])
    j.atmos[0].cloud[j.atmos[0].config.Cl['H2O']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['H2O']])
    j.atmos[0].cloud[j.atmos[0].config.Cl['NH4SH']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['NH4SH']])
    j.atmos[0].cloud[j.atmos[0].config.Cl['NH3']] = np.zeros_like(j.atmos[0].cloud[j.atmos[0].config.Cl['NH3']])

    rb.atm_modify.model_tceq(j, plotting=True)


    References
    ------------
    
    Todo
    ----- 

    Notes
    -------
    12/16/19, CM, initial comit 
    4/2/2021 CM, updated H2S reference value to be consistent with Solar abundances 
    '''

    P = self.atmos[0].gas[self.atmos[0].config.C['P']]
    T = self.atmos[0].gas[self.atmos[0].config.C['T']]

    # Import the difference between dry and wet adiabat 

    # if adiabat.lower() == 'dry':
    #     T = TPprofile(P,adiabat='dry', updated=True)
    # elif adiabat.lower() == 'wet':
    #     T = TPprofile(P,adiabat='wet', updated=True)
    # elif adiabat>0 and adiabat < 1.5: 

    #     Td  = TPprofile(adiabat='dry') 
    #     Tw  = TPprofile(adiabat='wet') 
    #     T = Tw + (Td-Tw)*adiabat 
    # else:
    #     try: 
    #         T = j.atmos[0].gas[j.atmos[0].config.C['T']]
    #     except TypeError: 
    #         T = j.atmos.gas[j.atmos.config.C['T']]

    
    nh3_prof = np.copy(self.atmos[0].gas[self.atmos[0].config.C['NH3']])
    h2s_prof = np.copy(self.atmos[0].gas[self.atmos[0].config.C['H2S']])
    h2o_prof = np.copy(self.atmos[0].gas[self.atmos[0].config.C['H2O']])

    # Adjust deep_nh3 for ammonia 
    deep_nh3 *=(1-nh3_h2o)*N2H
    deep_h2s *= S2H 
    # Deep NH3 is the well mixed ammonia at depth! For calculation mid and deep are the same
    mid_nh3 = (deep_nh3 - (np.log(Pb_vm)-np.log(Pt_vm))*dnh3dlP*N2H)
    
    # Compute H2S based on reference value and how much ammonia is depleted at the top of the atmosphere
    mid_h2s = deep_h2s * mid_nh3/deep_nh3


    # Temporary solution for water clouds 
    rh_nh4sh_nh3 = 1.0
    rh_nh4sh_h2s = 1.0

    # Set the stratosphere abudance to 1e-13 
    idx_stratosphere = np.where(P > 0.2)[0][0]
    nh3_prof[0:idx_stratosphere] = 1e-13 

    for i in range(idx_stratosphere,len(P)):

        # Calculate the saturation pressure for the relevant gases 
        # --------------------------------------------------------------
        # Ammonia 
        psat_nh3 = chemistry.ConstituentProperties('NH3').Psat(T[i])
        pp_nh3= P[i] * mid_nh3

        # H2S 
        psat_h2s = chemistry.ConstituentProperties('H2S').Psat(T[i])
        pp_h2s= P[i] * mid_h2s

        # NH4SH 
        psat_nh4sh = chemistry.ConstituentProperties('NH4SH').Psat(T[i])
        # Modify the saturation pressure by the relative humidity 
        psat_nh4sh /= rh_nh4sh

        # H2O profile 
        psat_h2O = chemistry.ConstituentProperties('H2O').Psat(T[i])
        pp_h2o = P[i] * mid_h2s

        # Calculate the fraction of depletion due to clouds (See Tollefson Thesis, Section 5.3) 
        # --------------------------------------------------------------
        a = 1 
        b = -(1+mid_nh3/mid_h2s)
        c =  mid_nh3/mid_h2s - psat_nh4sh/(mid_h2s*P[i])**2


        z = (-b - np.sqrt(b**2 - 4*a*c))/(2*a)

        pp_h2s = mid_h2s*P[i]
        pp_nh3 = mid_nh3*P[i]
        
        if z<0:
            z = 0 

        # Make the NH3 clouds, abundance is set by the abundance left after the NH4SH clouds 
        # --------------------------------------------------------------
        cpp_nh3 = (mid_nh3 -  z*mid_h2s)*P[i]

        # The relative humidity determines if over or undersaturate, and then changes the saturation pressure accordingly 
        if  cpp_nh3 > psat_nh3*rh_nh3:
            nh3_prof[i] = rh_nh3 * psat_nh3/P[i]
        else:
            # If we are ouyt
            nh3_prof[i] = mid_nh3
            h2s_prof[i] = deep_h2s

        #  NH4SH cloud, if partial pressure exceeds reaction pressure  
                # Make the NH3 clouds, abundance is set by the abundance left after the NH4SH clouds 
        # --------------------------------------------------------------
        if pp_h2s * pp_nh3 > psat_nh4sh and  cpp_nh3 < psat_nh3*rh_nh3:
            # print('We are making clouds at Pressure {:2.2f}'.format(P[i]))
            # print('nh3 ', nh3_prof[i],'z ',  z)
            nh3_prof[i] = rh_nh4sh_nh3 * (mid_nh3 - z*mid_h2s)
            h2s_prof[i] = rh_nh4sh_h2s * (mid_h2s - z*mid_h2s)


    # Mixing regime 
    idx_bottom  = np.where(P < Pb_vm)[0][-1]
    idx_top     = np.where(P > Pt_vm)[0][0]
    nh3_prof[idx_top:idx_bottom] = nh3_prof[idx_top] + dnh3dlP*(np.log(P[idx_top:idx_bottom])-np.log(P[idx_top]))*N2H
    nh3_prof[idx_bottom:] = deep_nh3

    # If it drops below 0, send a warning 
    if np.any(nh3_prof<0): 
        print("Warning - Abundance dropped below 0. Forcing 0")
        nh3_prof[nh3_prof<0] = 0 

    # --------------------------------------------------------------
    # Water abundance 

    # Set the stratosphere abudance to 1e-13 
    idx_cloudtop = np.where(P > 0.2)[0][0]
    h2o_prof[0:idx_cloudtop] = 1e-13 
    for i in range(len(P)): 
        # Caclulate if water should condense out  
        # H2O profile 
        psat_h2o = chemistry.ConstituentProperties('H2O').Psat(T[i])
        pp_h2o = P[i] * deep_h2o*O2H 
        # Make the water clouds 
        if  pp_h2o > psat_h2o*rh_h2o:
            h2o_prof[i] = rh_h2o * psat_h2o/P[i]
            P_h2o = P[i] # Save pressure at which the water vapor is condensing out 
        else: 
            h2o_prof[i] = deep_h2o*O2H 



    # --------------------------------------------------------------
    # make H2O cloud, only triggers once and introduces a step function  
    idx_h2o = np.argmax(P>P_h2o)
    nh3_prof[idx_h2o:] += nh3_prof[idx_h2o]*nh3_h2o 


    if plotting: 
        fig,ax = plt.subplots(figsize=(9,9)) 
        ax.plot(self.atmos[0].gas[self.atmos[0].config.C['NH3']]/N2H,P,color='gray', label=r'NH$_3$')
        ax.plot(self.atmos[0].gas[self.atmos[0].config.C['H2S']]/S2H,P,color='yellow', label=r'H$_2$S')
        ax.plot(self.atmos[0].gas[self.atmos[0].config.C['H2O']]/O2H,P,color='navy', label=r'H$_2$O')    


        ax.plot(nh3_prof/N2H,P,linestyle='--',color='gray', label=r'Modified NH$_3$')
        ax.plot(h2s_prof/S2H,P,linestyle='--',color='yellow', label=r'Modified H$_2$S')
        ax.plot(h2o_prof/O2H,P,linestyle='--',color='navy', label=r'Modified H$_2$O')

        ax.set_yscale('log')
        ax.invert_yaxis()
        ax.set_ylim([Pb_vm,0.1])

        ax.set_ylabel('P [bar]')
        ax.set_xlabel('Abundance [solar]')
        plt.legend()
        plt.show()

    # Updates the ammonia abundance 
    self.atmos[0].gas[self.atmos[0].config.C['NH3']] = np.copy(nh3_prof)
    self.atmos[0].gas[self.atmos[0].config.C['H2S']] = np.copy(h2s_prof)
    self.atmos[0].gas[self.atmos[0].config.C['T']]   = np.copy(T)
    self.atmos[0].gas[self.atmos[0].config.C['H2O']] = np.copy(h2o_prof)

    return None

