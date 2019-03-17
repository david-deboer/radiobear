from radiobear import chemistry
import matplotlib.pyplot as plt


def modify(gas, cloud, C, Cl):

    comment = "'Standard' Neptune atmospheric tweaking from Imke's code..."
    # Ad hoc "options"
    option = 'A'  # The "normal" NH4SH option, other is "B"
    depletech4 = True
    ch4 = chemistry.ConstituentProperties('CH4')
    print(ch4)
    nAtm = len(gas[C['P']])
    deep_CH4 = 0.04
    for i in range(nAtm):
        Plyr = gas[C['P']][i]
        Tlyr = gas[C['T']][i]
        # # ## Process CH4
        if depletech4:
            if option == 'A':

                if Plyr > 1.3:
                    gas[C['CH4']][i] = deep_CH4

    for i in range(nAtm):
        j = nAtm - i - 1
        Plyr = gas[C['P']][j]
        Tlyr = gas[C['T']][j]
        # # ## Process CH4
        if depletech4:
            if option == 'A':
                Psat_CH4 = ch4.Psat(Tlyr)
                if Plyr * gas[C['CH4']][j] > Psat_CH4:
                    gas[C['CH4']][j] = Psat_CH4 / Plyr
                    cloud[Cl['CH4']][j] = 1.0

        # ## Process CO
        # ## Process CO
        gas[C['CO']][i] = 0.0
        # ## Process CO13
        gas[C['CO13']][i] = 1.0E-2 * gas[C['CO']][i]

        # ## Process HCN, SOLN, PH3
        gas[C['HCN']][i] = 0.0
        gas[C['SOLN']][i] = 0.0
        gas[C['PH3']][i] = 0.0

    return comment, gas, cloud
