from radiobear import planet

j = planet.Planet('jupiter', plot_atm=False, verbose=False)
j.run('1:10:1', save_alpha='memory')

catch_data = []
for s in range(10):
    scale = (s + 1) / 5.0
    print("scale {}".format(scale))
    catch_data.append(j.run('1:10:1', scale=scale, get_alpha='memory'))

for cd in catch_data:
    print(cd)
