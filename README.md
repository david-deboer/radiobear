RadioBEAR
========

Radio BErkeley Atmospheric Radiative-transfer (RadioBEAR)

planetary atmosphere code to calculate the brightness temperature of planetary
atmospheres in the meter-to-millimeter wavelength range.

If used, please reference
1. de Pater et al 2019:   de Pater, I., R. J. Sault, M. H. Wong, L. N. Fletcher, D. DeBoer,
B. Butler,  2019. Jupiter's ammonia distribution derived from VLA maps
at 3--37 GHz. Icarus, 322, 168-191.
2. de Pater et al.(2014)] de Pater, I.,  Fletcher, L. N., Luszcz-Cook, S., DeBoer, D., Butler, B., Hammel, H. B., Sitko, M. L., Orton, G., Marcus, P. S.}, Neptune's global circulation deduced from multi-wavelength observations. Icarus 237, 211

Installing:
1. Pull radioBEAR from github
2. From the top-level 'radiobear' directory, install (`pip install .`)*
3. Create/move to your working directory (so, outside of the radiobear installation)
4. Initialize your working directory with various directories/planet files
      - from within your working directory type `initial_planet_setup.py`
5. If you have other .par, tweak, atmosphere, etc files, move them to the planet sub-directory in the new working area

[*NB: python setup.py install does not work.]

Before you start edit the config file to be what you want:
1. Nearly all of the parameters are set within this configuration file
2. Default name is `config.par` in each of the planet sub-directories
3. You may use different files, just call planet with `config='filename'`, which must reside within <planet_name>
4. The defaults are set in `default_config.json`, which also sets up which parameters are contained within the config

Note that the defaults are in default_config.json and default_state.json.  The two files roughly break out different
"types", but work identically.

Note that you may also set any config parameter at runtime by setting `par=val`


Here is an example, from working area:
```
~/rbwork> from radiobear import planet
~/rbwork> j = planet.Planet('jupiter')
~/rbwork> catch_data = j.run(freqs='1:100:5', b='disc')
```

time-stamped data file is written to Output

run log file is written to Logs

'plot' is a keyword (see below)

`catch_data` is an instance of `DataReturn`, which holds the data you probably want.

options for name:  `Jupiter`, `Saturn`, `Uranus`, `Neptune`

examples for freqs:
* `freqs = 1.42`    ==> single frequency at 1.42 GHz
* `freqs = [1.4,2.5,5.5,8.4,12.1]`  ==> uses these frequencies (can be a list, csv string or numpy array)
* `freqs = '1:10:1'` ==> this will generate a range as start:stop:step (stop is always included)
* `freqs = '1;100;20'` ==> this will generate a log range as start;stop;nvalues
* `freqs = 'freq.dat'`   ==> (string within '') reads in those frequencies, one per line


examples for b:
* `b = 0.1`  ==> generates a full image at that resolution (see blocks)
* `b = 'stamp'` ==> generates a small image (queries for extents)
* `b = [[0.0,0.0],[0.1,0.0],...]`  ==> generates at listed points
* `b = [0.0,0.0]` ==> same as above at that one point
* `b = 'disc'` ==> disc-averaged brightness temperature
* `b = '0.0,0.2,0.4,0.6,0.8,1.0<45'` ==> will use those values as an angle of 45deg, (default angle is 0.0)
* `b = '0.0:1.0:0.1<45'` ==> range start:stop:step<angle                                                  >

========
 Routines for viewing the data files in the Output directory are in the 'plotting.plt' module.  

 To view spectra try:
```
  > from radiobear.plotting import plt
  > d=plt.Tb(legend=True)
```
See options under plt.Tb()
