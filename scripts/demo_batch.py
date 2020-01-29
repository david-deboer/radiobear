from radiobear import planet
import numpy as np

j = planet.Planet('jupiter', plot_atm=False, verbose=False)
cd = j.run('1:10:1', save_alpha='memory')  # could be save_alpha='file'
print("Calculate first")
print(cd)

nlayers = np.shape(j.atmos[0].gas)[1]

scale_ver = {}
scale_ver['const'] = 0.5
scale_ver['total'] = [4.0] * nlayers
scale_ver['NH3'] = {'nh3': [3.0] * nlayers}

for k, scale in scale_ver.items():
    print("Scale Version {}".format(k))
    cd = j.run('1:10:1', scale=scale, get_alpha='memory')  # could get get_alpha='file'
    print(cd)
