from radiobear import planet


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
scale_ver['total'] = [4.0] * nlayers
scale_ver['NH3'] = {'nh3': [3.0] * nlayers}
scale_ver['calc'] = make_random_scale(j.atmos[0], j.config.C)

for k, scale in scale_ver.items():
    print("Scale Version {}".format(k))
    cd = j.run('1:10:1', scale=scale, get_alpha='memory')  # could get get_alpha='file', could re save_alpha  # noqa
    print(cd)
