from radiobear import planet
import matplotlib.pyplot as plt


def make_random_scale(atm, cols):
    scale = {'nh3': [], 'h2o': []}
    for p in atm.gas[cols['P']]:
        if p > 100.0:
            scale['nh3'].append(1.0)
            scale['h2o'].append(2.0)
        elif p > 10.0:
            scale['nh3'].append((p - 5.0) / 50.0)
            scale['h2o'].append(0.5)
        elif p > 1.0:
            scale['nh3'].append(0.25)
            scale['h2o'].append((p - 0.9) / 8.0)
        else:
            scale['nh3'].append(2.0)
            scale['h2o'].append(0.0001)
    return scale


j = planet.Planet('jupiter', plot_atm=False, verbose=False)
cd = j.run('1:10:1', save_alpha='memory')  # could be save_alpha='file'
print("Calculate first")
print(cd)

nlayers = len(j.atmos[0].gas[0])

scale_ver = {}
scale_ver['const'] = 0.5
scale_ver['ramp_on_total'] = []
Pressures = j.atmos[0].gas[2]
for i in range(nlayers):
    p = Pressures[i]
    scale_ver['ramp_on_total'].append(2.0 * (Pressures[-1] - p / 1.5) / Pressures[-1])
scale_ver['calc_dict'] = make_random_scale(j.atmos[0], j.config.C)

for k, scale in scale_ver.items():
    print("Scale Version {}".format(k))
    cd = j.run('1:10:1', scale=scale, get_alpha='memory')  # could get get_alpha='file', could re save_alpha  # noqa
    print(cd)

plt.figure('Scales')
plt.semilogx([Pressures[0], Pressures[-1]], [scale_ver['const'], scale_ver['const']], label='const')
plt.semilogx(Pressures, scale_ver['ramp_on_total'], label='ramp_on_total')
plt.semilogx(Pressures, scale_ver['calc_dict']['nh3'], label='nh3')
plt.semilogx(Pressures, scale_ver['calc_dict']['h2o'], label='h2o')
plt.legend()
