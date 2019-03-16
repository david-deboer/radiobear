import sys
import os
import emcee
from shutil import copyfile
import numpy as np
import time
import EmceeFiles.emcee_input as emcee_input
import planet

### Load in real data and model spectrum

def get_obs_spectrum(specdatfile):
    '''Read spectrum data file
        from ./Output/.dat
        Should have format:
        wavelengths (GHz), Tb (K), Tb_err(K)'''

    spec = []
    with open(specdatfile, 'r') as tmp:
        for line in tmp:
            if not line.startswith("#"):
                values=[float(x) for x in line.split()]
                if len(spec) == 0:
                    spec = values
                else:
                    spec = np.vstack((spec, values))

    return spec

def get_model_spectrum(val, name, tweakmodule, tweakTmp, freqs, b, par_names, limits):
    '''Generate model spectrum'''

    update_tweak(tweakmodule,tweakTmp,par_names, val, limits)
    plan = planet.planet(name, plot=False)
    datFile = plan.run(freqs=freqs, b=b)

    spec = []
    with open(datFile, 'r') as tmp:
        for line in tmp:
            if not line.startswith("#"):
                values=[float(x) for x in line.split()]
                if len(spec) == 0:
                    spec = values
                else:
                    spec = np.vstack((spec, values))

    os.remove(datFile)
    return spec

def update_tweak(tweakmodule, tweakTmp, par_names, guess, limits):
    '''Update tweak file val with new val'''

    copyfile(tweakTmp, tweakmodule)
    for i in range(len(par_names)):
        f = open(tweakmodule, 'r')
        filedata = f.read()
        f.close()
        newdata = filedata.replace(par_names[i], str(guess[i]))
        f = open(tweakmodule, 'w')
        f.write(newdata)
        f.close()
    f.close()


### Compute probabilities

def lnprior(theta,limits): #flat priors
    tmp = 1
    for i in range(len(theta)):
        if limits[i][0] < theta[i] < limits[i][1]:
            tmp*=1
        else: tmp*=0
    if tmp==1: return 0.0
    else: return -np.inf

def lnlike(theta, x, y, yerr, name, tweakmodule, tweakTmp, freqs, b, par_names, limits):
    parvals=theta
    spec = get_model_spectrum(parvals, name, tweakmodule, tweakTmp, freqs, b, par_names, limits)
    ymodel = spec[:,1]

    sigsq=yerr**2
    lnP=-0.5*np.sum((y-ymodel)**2./sigsq+np.log(2*np.pi*sigsq))
    return lnP


def lnprob(theta,x,y,yerr,name,tweakmodule,tweakTmp,freqs,b,par_names, limits):
    lp = lnprior(theta,limits)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,x,y,yerr,name,tweakmodule,tweakTmp,freqs,b,par_names, limits)

def run_emcee_spectrum(sampler,pos, nsteps,outdatfile,lnprobfile=None):
    t0=time.time()
    for result in sampler.sample(pos, iterations=nsteps, storechain=True):
        position = result[0]
        f = open(outdatfile, "a")
        #pdb.set_trace()
        for k in range(position.shape[0]):
            f.write("{0:4d} {1:s}\n".format(k, " ".join(map(str,position[k]))))
        f.close()
        if lnprobfile is not None:
            probability = result[1]
            f = open(lnprobfile,"a")
            for k in range(probability.shape[0]):
                f.write("{0:4d} {1:s}\n".format(k, repr(probability[k])))
            f.close()
    t1=time.time()
    print ('run took {0} minutes to travel {1} steps'.format((t1-t0)/60.,nsteps))
    print ('acceptance fraction:',sampler.acceptance_fraction)

    return sampler

def run_emcee_spectrum_new():
    ##name='Neptune', refData='./Neptune/refSpectrum.dat', tweakmodule = './Neptune/NeptuneTweak_mcmc.py', tweakTmp = './Neptune/tmp.py'
    #guess=1.0, limits=[0.75,1.25],
    #    nwalker=10, threads=1, nsteps=100,
    #        outdatfile='out.dat',lnprobfile='lnprobfile.dat',lnprobfn=lnprob

    inp = emcee_input.gen_emcee_input()

    name = inp['planet']
    refData = inp['refData']
    tweakmodule = inp['tweakmodule']
    tweakTmp = inp['tweakTmp']

    freqs = inp['freqs']
    b = inp['b']

    par_names = inp['parameters']['names']
    guess = inp['parameters']['guesses']
    limits = inp['parameters']['limits']
    nwalker = inp['nwalkers']
    threads = inp['threads']
    nsteps = inp['nsteps']

    outdatfile = inp['outdatafile']
    lnprobfile = inp['lnprobfile']

    #Check if we're going to overwrite a file
    if os.path.isfile(outdatfile):
        overwrite=False
        if overwrite is False:
            print ("Please choose another filename and rerun: {}".format(outdatfile))
            return
    continueprob, newprob= check_likefile(lnprobfile, "lnprobfile")
    if continueprob is False:
        print ("Exiting because of a problem with lnprobfile.")
        return
    if newprob is False:
        print ("lnprobfile already exists! Exiting")
        return
    f = open(lnprobfile,"w")
    f.close()

    name = str.capitalize(name)
    ndim = len(guess)

    obs_spec =  get_obs_spectrum(refData)
    x = obs_spec[:,0]
    y = obs_spec[:,1]
    yerr = obs_spec[:,2]

    pos = [guess + 1e-4*np.random.randn(ndim) for i in range(nwalker)]

    sampler = emcee.EnsembleSampler(nwalker, ndim, lnprob, args=(x,y,yerr,name,tweakmodule,tweakTmp,freqs,b, par_names,limits), threads=threads)

    f = open(outdatfile, "w")
    f.close()
    sampler=run_emcee_spectrum(sampler, pos, nsteps, outdatfile, lnprobfile=lnprobfile)
    return sampler

def run_emcee_spectrum_append():

    inp = emcee_input.gen_emcee_input()

    name = inp['planet']
    refData = inp['refData']
    tweakmodule = inp['tweakmodule']
    tweakTmp = inp['tweakTmp']

    par_names = inp['parameters']['names']
    guess = inp['parameters']['guesses']
    limits = inp['parameters']['limits']
    nwalker = inp['nwalkers']
    threads = inp['threads']
    nsteps = inp['nsteps']

    freqs = inp['freqs']
    b = inp['b']


    datfile = inp['outdatafile']
    lnprobfile = inp['lnprobfile']

    name = str.capitalize(name)


    obs_spec =  get_obs_spectrum(refData)
    x = obs_spec[:,0]
    y = obs_spec[:,1]
    yerr = obs_spec[:,2]

    if os.path.isfile(datfile) is False:
        print ("{0} is not found. Did you mean to run run_emcee_spectrum_new?".format(datfile))
        return
    append = True
    print ("appending to file {0}?".format(datfile))
    if append is False:
        print ("If you want to continue running emcee on {0}, then say so! Returning.".format(datfile))
        return
    # Check probability file
    continueprob, newprob= check_likefile(lnprobfile, "lnprobfile")
    if continueprob is False:
        print ("Exiting because of a problem with lnprobfile.")
        return
    if newprob is True:
        print ("Warning: printing probabilities to a NEW file, while appending walker positions to an OLD file.")
        f = open(lnprobfile,"w")
        f.close()

    samples,ndim,nwalkers=read_emcee_datfile(datfile)
    pos=samples[-1,:,:]
    sampler = emcee.EnsembleSampler(nwalker, ndim, lnprob, args=(x,y,yerr,name,tweakmodule,tweakTmp,freqs,b, par_names,limits), threads=threads)
    sampler=run_emcee_spectrum(sampler, pos, nsteps, datfile,lnprobfile=lnprobfile)
    return sampler


def check_likefile(filename, filestring):
    new=False
    if filename is None:
        continue_run=False
    else:
        if os.path.isfile(filename) is False:
            print ('writing new file')
            continue_run=True
            new= True
        else:
            continue_run=True
            print ('appending to existing file')
    return continue_run, new


def read_emcee_datfile(datfile):
    out=[]
    with open(datfile,'r') as tmp:
        for line in tmp:
            values=[float(x) for x in line.translate('[]').split()]
            if len(out) == 0:
                out=values
            else:
                out=np.vstack((out,values))

    nwalkers=np.max(out[:,0]) +1
    ndim=out.shape[1]-1
    samples=out[:,1:].reshape(-1,nwalkers,ndim)
    return samples,ndim,nwalkers

def read_emcee_probfile(probfile):
    file=open(probfile,'r')
    tmp=file.readlines()
    file.close()
    out=[]
    for line in tmp:
        values=[float(x) for x in line.split()]
        if len(out) == 0:
            out=values
        else:
            out=np.vstack((out,values))

    nwalkers=np.max(out[:,0]) +1
    probabilities=out[:,1].reshape(-1,nwalkers)
    return probabilities,nwalkers
