import numpy as np
import emcee
import scipy.optimize as op
import corner

def load_model(datfile):
    out = []
    with open(datfile, 'r') as tmp:
        for line in tmp:
            if not line.startswith("#"):
                values = [float(x) for x in line.translate('[]').split()]
                if len(out) == 0:
                    out=values
                else:
                    out = np.vstack((out,values))
    return out

def lnlike(theta, x, y, yerr):
    '''Generate likelihood estimation.
        Gaussian underestimated by a fractional amount lnf'''
    m, b, lnf = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    return -0.5*(np.sum((y - model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
    '''Generate parameter prior.
        Uniform distribution.'''
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf

def lnprob(theta, x, y, yerr):
    '''Full probability function'''
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

## Choose the "true" parameters.
#m_true = -0.9594
#b_true = 4.294
#f_true = 0.534

# Generate some synthetic data from the model.
#N = 50
#x = np.sort(10*np.random.rand(N))
#yerr = 0.1+0.5*np.random.rand(N)
#y = m_true*x+b_true
#y += np.abs(f_true*y) * np.random.randn(N)
#y += yerr * np.random.randn(N)

## Initialize guess (result), which is based on typical least-squares optimizaiton
#nll = lambda *args: -lnlike(*args)
#result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
#m_ml, b_ml, lnf_ml = result["x"]

## Set up emcee with walkers and their inital positions
#ndim, nwalker = 3, 100
#pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalker)]
#sampler = emcee.EnsembleSampler(nwalker, ndim, lnprob, args=(x,y,yerr))

## Run code and plot/print results
#sampler.run_mcmc(pos, 500)
#samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

#fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
#                    truths=[m_true, b_true, np.log(f_true)])
#fig.savefig("triangle.png")

#samples[:, 2] = np.exp(samples[:, 2])
#m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
#                             zip(*np.percentile(samples, [16, 50, 84],
#                                                axis=0)))

#print (m_mcmc, b_mcmc, f_mcmc)
