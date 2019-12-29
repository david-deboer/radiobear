#! /usr/bin/env python

from __future__ import print_function
import os
import shutil
from radiobear import utils as rbu

path_to_radiobear_planets = rbu.get_location_for_radiobear_setup()

planets_to_update = ['Jupiter', 'Saturn', 'Uranus', 'Neptune']
files_to_update = ['bmap.py']

for planet in planets_to_update:
    radiobear_path = os.path.join(path_to_radiobear_planets, planet)
    for pf in files_to_update:
        print("Moving {} to {}".format(pf, planet))
        init_location = os.path.join(radiobear_path, pf)
        final_location = os.path.join(planet, pf)
        shutil.copy(init_location, final_location)
