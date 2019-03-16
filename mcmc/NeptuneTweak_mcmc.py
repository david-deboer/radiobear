def modify(gas,cloud,C,Cl):

    comment = "Tweak Neptune CH4 content by amount given at all layers by mcmc code - for testing purposes"

    nAtm = len(gas[C['P']])
    for i in range(nAtm):
        Plyr = gas[C['P']][i]
        Tlyr = gas[C['T']][i]

        ### Process CH4
 
        ### Process H2S

        if Plyr < 200 and Plyr > 50:
            gas[C['H2S']][i] = 0.661231765056*gas[C['H2S']][i]
        elif Plyr < 50:
            gas[C['H2S']][i] = 0.603937377839*gas[C['H2S']][i]
        
        ### Process NH3


        ### Process CO
        gas[C['CO']][i] = 0.0
        ### Process CO13
        gas[C['CO13']][i] = 1.0E-2*gas[C['CO']][i]
        ### Process HCN, SOLN, PH3
        gas[C['HCN']][i] = 0.0
        gas[C['SOLN']][i] = 0.0
        gas[C['PH3']][i] = 0.0

    return comment, gas, cloud

