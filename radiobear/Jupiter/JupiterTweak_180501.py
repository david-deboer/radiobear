def modify(gas, cloud, C, Cl):

    comment = "Jupiter tweaking"

    nAtm = len(gas[C['P']])
    for i in range(nAtm):
        Plyr = gas[C['P']][i]
        Tlyr = gas[C['T']][i]

        if Plyr < 1e-3:
            gas[C['NH3']][i] = 1e-12
        # ## Process H2S

        # ## Process NH3

        # ## Process CO
        gas[C['CO']][i] = 0.0

        # ## Process CO13
        gas[C['CO13']][i] = 1.0E-2 * gas[C['CO']][i]

        # ## Process HCN, PH3
        gas[C['HCN']][i] = 0.0
        gas[C['PH3']][i] = 0.0

    return comment, gas, cloud
