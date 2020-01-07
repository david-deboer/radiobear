from radiobear import planet
from numpy import np

j = planet.Planet('jupiter', plot_atm=False, verbose=False)
j.run('1:10:1', save_alpha='memory')
nlayers = np.shape(j.atmos[0].gas)[1]

scale_ver = {}

scale_ver['const'] = (range(10) + 1) / 5.0
scale_ver['ramp'] = ['array of nlayers', 'another array of nlayers']
scale_ver['NH3'] = [{'NH3': 'array of nlayers'}, {'NH3': 'another array of nlayers'}]

use_scale = 'const'
catch_data = []
print('Using test scale {}'.format(use_scale))
for i, scale in enumerate(scale_ver[use_scale]):
    print("Version {}".format(i))
    print("scale {}".format(scale))
    catch_data.append(j.run('1:10:1', scale=scale, get_alpha='memory'))

for cd in catch_data:
    print(cd)
