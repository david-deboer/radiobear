from radiobear import planet
import numpy as np

j = planet.Planet('jupiter', plot_atm=False, verbose=False)
print("Computing first atmosphere to save into a file.")
j.run('1:10:1', save_alpha='file')
print("Now get into memory (you don't need to do file first - this is just for demo)")
j.run('1:10:1', get_alpha='file', save_alpha='memory')

nlayers = np.shape(j.atmos[0].gas)[1]

scale_ver = {}

scale_ver['const'] = np.arange(0.5, 2.1, 0.5)
scale_ver['ramp'] = ['array of nlayers', 'another array of nlayers']
scale_ver['NH3'] = [{'NH3': 'array of nlayers'}, {'NH3': 'another array of nlayers'}]

use_scale = 'const'
print('Using test scale {}'.format(use_scale))
for i, scale in enumerate(scale_ver[use_scale]):
    print("Scale Version {}".format(i))
    print(scale, type(scale))
    cd = j.run('1:10:1', scale=scale, get_alpha='memory')
    print(cd)
