from __future__ import print_function

import chemistry


def modify(gas, cloud, C, Cl):

    comment = "Jupiter tweaking - add in PH3"

    ph3 = chemistry.ConstituentProperties('PH3')
    deep_PH3 = ph3.solar

    nAtm = len(gas[C['P']])
    for j in range(nAtm):
        i = nAtm - j - 1
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
        gas[Cl['PH3']][i] = 0.0
        if j == 0:
            gas[C['PH3']][i] = deep_PH3
        else:
            gas[C['PH3']][i] = gas[C['PH3']][i + 1]
        Psat_PH3 = ph3.Psat(Tlyr)
        if Plyr * gas[C['PH3']][i] > Psat_PH3:
            gas[C['PH3']][i] = Psat_PH3 / Plyr
            gas[Cl['PH3']][i] = 1.0

    return comment, gas, cloud
