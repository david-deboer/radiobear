def modify(gas, cloud, C, Cl):

    comment = "No Tweaking done, except add in 0.0's"

    nAtm = len(gas[C['P']])
    for i in range(nAtm):

        # ## Process CO
        gas[C['CO']][i] = 0.0
        # ## Process CO13
        gas[C['CO13']][i] = 1.0E-2 * gas[C['CO']][i]

        # ## Process HCN, SOLN, PH3
        gas[C['HCN']][i] = 0.0
        gas[C['SOLN']][i] = 0.0
        gas[C['PH3']][i] = 0.0

    return comment, gas, cloud
