from __future__ import print_function, absolute_import, division

import chemistry


def get_value(p, t, q, c, v):
    """
    Put all the smart stuff here.  Can it be parameterized by a few good MCMC-able things?
    This is a sample from Josh...
    """

    gas = chemistry.ConstituentProperties(c)
    Psat_gas = gas.Psat(t)

    if c.upper() == 'H2S':
        if p < 43. and p * q * v > Psat_gas:  # Pressure greater than saturation pressure
            return str(1.0)
        elif p < 43. and p * q * v < Psat_gas:
            return str(v)
        else:
            return str(0.8)
    else:
        return str(1.0)
