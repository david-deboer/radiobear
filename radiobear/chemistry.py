import math
GravConst = 6.6738e-20  # km3/kg/s2
R = 8.314462            # Universal gas constant [J/K/mol]
AMU = 1.66056E-24       # [g]


def NULL(T=None, f=None):
    return 0.0


# atomic mass units [g/mol]
amu = {'H2': 2.016, 'HE': 4.003, 'H2S': 34.076, 'NH3': 17.030, 'H2O': 18.015,
       'CH4': 16.04, 'PH3': 33.997, 'NH4SH': 51.110, 'NULL': 0.0}

# solar abundances (taken from dePater/Lissauer p 80)
solar_abundance = {'H2': 0.835, 'HE': 0.195, 'H2O': 1.70E-3, 'CH4': 7.94E-4, 'NH3': 2.24E-4,
                   'H2S': 3.70E-5, 'PH3': 7.50E-7, 'NULL': 0.0}

# triple points [K]
triple_point = {'H2S': 187.61, 'NH3': 195.5, 'CH4': 90.7, 'PH3': 0.0, 'H2O': 273.16, 'NULL': 0.0}

# specific heat at constant pressure over R.  High-T limit for H2
specific_heat = {'H2': 3.5, 'HE': 2.50, 'H2O': 4.00, 'NH3': 4.46, 'H2S': 4.01, 'CH4': 4.50, 'NULL': 0.0}


def REFR_H2O(T, f=1.0):
    if T is None:
        T = 100.0
    return (245.0 + 1.28E6 / T)  # Janssen p218--or is the h2o cloud really almost an ocean,see liquid water p298 UFM*/


# refractivity coefficients (Rc) R = Rc*P*(293/T) then n = (R/1E6) + 1. Average values of microwave
refractivity = {'H2': 124.43, 'HE': 35.832, 'H2O': REFR_H2O, 'H2S': 2247.0, 'NH3': 2700.0, 'CH4': 413.0, 'NULL': NULL}


def NH3_over_NH3_ice(T):  # Briggs and Sackett
    a1 = -4122.0
    a2 = 27.8632
    a3 = -1.8163
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def NH3_over_liquid_NH3(T):  # Briggs and Sackett
    a1 = -4409.3512
    a2 = 63.0487
    a3 = -8.4598
    a4 = 5.51E-3
    a5 = 6.80E-6
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def H2S_over_H2S_ice(T):  # Allen, Giauque/Blue
    a1 = -2920.6
    a2 = 14.156
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def H2S_over_liquid_H2S(T):  # Allen, Giauque/Blue
    a1 = -2434.62
    a2 = 11.4718
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def H2O_over_ice(T):  # Briggs and Sackett
    a1 = -5631.1206
    a2 = -22.1791
    a3 = 8.2312
    a4 = -3.861449e-2
    a5 = 2.77494e-5
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def H2O_over_water(T):  # Briggs and Sackett
    a1 = -2313.0338
    a2 = -177.848
    a3 = 38.053682
    a4 = -0.13844344
    a5 = 7.4465367e-5
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def CH4_over_CH4_ice(T):  # dePater and Massie
    a1 = -1168.1
    a2 = 10.710
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def CH4_over_liquid_CH4(T):  # dePater and Massie
    a1 = -1032.5
    a2 = 9.216
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def NH4SH(T):  # Lewis
    a1 = -10834.0
    a2 = 34.151
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


def PH3_over_PH3_ice(T):  # Orton/Kaminski
    a1 = -1830.0
    a2 = 9.8225
    a3 = 0.0
    a4 = 0.0
    a5 = 0.0
    sp = a1 / T + a2 + a3 * math.log(T) + a4 * T + a5 * T * T
    return math.exp(sp)


Psat_all = {'H2S': {'ice': H2S_over_H2S_ice, 'liquid': H2S_over_liquid_H2S},
            'NH3': {'ice': NH3_over_NH3_ice, 'liquid': NH3_over_liquid_NH3},
            'CH4': {'ice': CH4_over_CH4_ice, 'liquid': CH4_over_liquid_CH4},
            'H2O': {'ice': H2O_over_ice, 'liquid': H2O_over_water},
            'PH3': {'all': PH3_over_PH3_ice},
            'NH4SH': {'all': NH4SH},
            'NULL': {'all': NULL}}


class ConstituentProperties:
    def __init__(self, constituent):
        self.constituent = constituent.upper()
        self.chkey = {'amu': 'NULL', 'triple_point': 'NULL', 'cp': 'NULL', 'solar': 'NULL',
                      'psat': 'NULL', 'refractivity': 'NULL'}

        if self.constituent in amu:
            self.chkey['amu'] = self.constituent
        self.amu = amu[self.chkey['amu']]

        if self.constituent in triple_point:
            self.chkey['triple_point'] = self.constituent
        self.triple_point = triple_point[self.chkey['triple_point']]

        if self.constituent in specific_heat:
            self.chkey['cp'] = self.constituent
        self.cp = specific_heat[self.chkey['cp']]

        if self.constituent in solar_abundance:
            self.chkey['solar'] = self.constituent
        self.solar = solar_abundance[self.chkey['solar']]

        if self.constituent in Psat_all:
            self.chkey['psat'] = self.constituent
        self.Psat_func = Psat_all[self.chkey['psat']]

        if self.constituent in refractivity:
            self.chkey['refractivity'] = self.constituent
        self.refractivity_func = refractivity[self.chkey['refractivity']]

    def Psat(self, T):
        k = 'all'
        if len(self.Psat_func.keys()) > 1:
            if T < self.triple_point:
                k = 'ice'
            else:
                k = 'liquid'
        return self.Psat_func[k](T)

    def refractivity(self, T=None, f=None):
        if type(self.refractivity_func) == float:
            return self.refractivity_func
        return self.refractivity_func(T=T, f=f)
