RadioBEAR
========

Radio BErkeley Atmospheric Radiative-transfer (RadioBEAR)

planetary atmosphere code to calculate the brightness temperature of planetary
atmospheres in the meter-millimeter-wave regions.

KNOWN ISSUES:
- For some matplotlib installs, the `plot_bright` plots sometimes cause the terminal to freeze.  

Installing:
1. Download radioBEAR from github
2. From the top-level 'radiobear' directory, install (`pip install .`)
3. Create/move to your working directory (so, outside of the radiobear installation)
4. Initialize your working directory with various directories/planet files
      - from within your working directory type `initial_planet_setup.py`
5. If you have other .par, tweak, atmosphere, etc files, move them from where they were into the new planet folder.


Before you start:
1. Set up your **config file**
    1. Nearly all of the parameters are set within this configuration file
    2. The defaults are set in `default_config.json`, which also sets up which parameters are contained within the config
    3. Default name is `config.par` in each of the planet sub-directories
    4. You may use different files, just call planet with `config='filename'`, which must reside within <planet_name>
    5. The tweakmodule filename gets set within the config, which adjusts the read in atmosphere - make sure it is what you want or make it so.
    6. For mcmc, the scalemodule also gets set here (and the scalefilename)
2. When you make a planet instance you may also set state_variables, they are initialized in state_variables.py
3. When you make a planet instance you may specify a mode, which sets the state_variables to various configurations (see state_variables.py)
4. The planet call variables and defaults are shown below.  kwargs may be one of the state_variables
5. By convention, the instance is set to the lower-case first letter of the planet (not required)


Within a python environment here is an example:
```
> import radiobear as rb
> j = rb.planet.Planet('jupiter', plot=False)
> catch_data = j.run(freqs='1:100:5', b='disc')
```
time-stamped data file is written to Output and log file to Logs
'plot' is a keyword (see below)

The declaration for planet is
```
class Planet:
    def __init__(self, name, mode='normal', config='config.par', **kwargs)
```

options for name:  Jupiter, Saturn, Uranus, Neptune

options for mode:  normal, mcmc, batch, use_alpha, scale_alpha

options for freqs:
* `freqs = 1.42`    ==> single frequency at 1.42 GHz
* `freqs = [1.4,2.5,5.5,8.4,12.1]`  ==> uses these frequencies (can be a list, csv string or numpy array)
* `freqs = '1:10:1'` ==> this will generate a range as start:stop:step (stop is always included)
* `freqs = '1;100;20'` ==> this will generate a log range as start;stop;nvalues
* `freqs = 'freq.dat'`   ==> (string within '') reads in those frequencies, one per line


options for b:
* `b = 0.1`  ==> generates a full image at that resolution (see blocks)
* `b = 'stamp'` ==> generates a small image (queries for extents)
* `b = [[0.0,0.0],[0.1,0.0],...]`  ==> generates at listed points
* `b = [0.0,0.0]` ==> same as above at that one point
* `b = 'disc'` ==> disc-averaged brightness temperature
* `b = '0.0,0.2,0.4,0.6,0.8,1.0<45'` ==> will use those values as an angle of 45deg, (default angle is 0.0)
* `b = '0.0:1.0:0.1<45'` ==> range start:stop:step<angle                                                  >

output of <>.state() (kwargs)

radiobear.planet.Planet state variables and defaults
	batch_mode:  False           # run in batch_mode
	initialize:  True            # initialize run when class instantiated
	write_output_files:  True    #
	write_log_file:  True.       #
	plot:  None                  # will set both plot_atm and plot_bright, only used when class set up
	plot_atm:  False             # plot the atmosphere data
	plot_bright:  True           # plot the radiative transfer data
	verbose:  True               #
	generate_alpha:  False       # generate an absorption profile file
	use_existing_alpha:  False.  # reuse the generated file
	scale_existing_alpha:  False # scale and use the generated file
	normalize_weighting:  True.  # normalize to 1
	output_type:  frequency.     # frequency or wavelength
	log_directory:  Logs.        # sub-directory for logs
	output_directory:  Output.   # sub-directory for output
	scratch_directory:  Scratch. # sub-directory for eg. Absorption etc.
