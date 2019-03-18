#! /usr/bin/env python
from __future__ import print_function
import os
import shutil
from radiobear import utils as rbu

path_to_planets = rbu.get_location_for_planet_setup()

planets_to_setup = ['Jupiter', 'Saturn', 'Uranus', 'Neptune']

for planet in planets_to_setup:
    print("Setting up {}".format(planet))
    if not os.path.isdir(planet):
        os.mkdir(planet)
    else:
        print("{} directory exists, overwriting within it.".format(planet))
    planet_path = os.path.join(path_to_planets, planet)
    for pf in os.listdir(planet_path):
        if pf[0] in ['.', '_']:
            continue
        init_location = os.path.join(planet_path, pf)
        final_location = os.path.join(planet, pf)
        shutil.copy(init_location, final_location)

other_directories = ['Logs', 'Output', 'Scratch']
for other in other_directories:
    if not os.path.isdir(other):
        print("Making sub-directory {}".format(other))
        os.mkdir(other)
