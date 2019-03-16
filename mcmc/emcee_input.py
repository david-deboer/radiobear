import numpy as np


def gen_emcee_input():
    planet = 'Neptune'
    # refData - reference spectrum; should have three columns: frequencies, Tb, and Tb_error
    # scaleFile - name of tweak file read by config.par
    refData = './emceeTest/refSpectrum.dat'

    b = 'disc'

    # ## Emcee Inputs
    # Parameters:
    # names - list of identifying names for the free parameters that emcee guesses and replaces in the scale file. Equal to fractional H2S depletion at two pressure ranges in this example. Names should be ALL CAPS.
    # guesses = list of initial guesses to the free parameters listed in 'names' to change in scale file - equal to the equilibrium value of H2S in this example
    # limits = list of values walkers are limited to search - between 0% to 120% abundance of H2S in this example
    # scalefile = name of the initial scale file for those constituents.  If it doesn't exist it will be initialized to 'guesses'
    # nwalker = number of Goodman & Weare walkers
    # threads = if > 1 will implement parallelization
    # nsteps = number of iterations to take before stopping - may append to prior runs by running emcee_append
    parameters = {'names':    ['H2S'],
                  'guesses':  [0.9],
                  'limits':   [[0.25, 1.75]]}

    nwalkers = 10
    threads = 1
    nsteps = 100

    # ## Data output files for coninuing or analyzing emcee runs. Files must not exist prior to running a new emcee run.
    outdatafile = './emceeTest/out.dat'
    lnprobfile = './emceeTest/lnprobfile.dat'
    configfile = 'config.par'

    # extract frequencies from the data
    refd = np.genfromtxt(refData, comments='#')
    freqs = sorted(list(refd.T[0]))

    # Make everything a dict
    inp = {
           'planet':        planet,
           'refData':       refData,
           'freqs':         freqs,
           'b':             b,
           'parameters':    parameters,
           'nwalkers':      nwalkers,
           'threads':       threads,
           'nsteps':        nsteps,
           'outdatafile':   outdatafile,
           'lnprobfile':    lnprobfile,
           'configfile':    configfile}

    print(inp['parameters'])
    return inp
