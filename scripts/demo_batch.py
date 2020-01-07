from radiobear import planet

j = planet.Planet('jupiter', plot_atm=False)
j.run('1:10:1', save_alpha='memory')

for s in range(10):
    scale = s / 5.0
    j.run('1:10:1', scale=scale, get_alpha='memory')
