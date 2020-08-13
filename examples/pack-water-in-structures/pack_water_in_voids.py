""" In this example, we first fetch a block of amorphous silica of size (357Å, 143Å, 143Å). We then carve spheres of random sizes at random locations, and finally pack 10,000 water molecules into the voids. 
"""

import numpy import random
from molecular_builder import fetch_prepared_system, carve_geometry, pack_water
from molecular_builder.geometry import SphereGeometry, BoxGeometry

# Amorphous silica of size (357, 143, 143)
amorph = fetch_prepared_system("amorphous_silica_1")

num_spheres = 20
for sphere in range(num_spheres):
    i, j, k, l = random.uniform(size=4)
    x, y, z, r = i*357, j*143, k*143, l*50
    geometry = SphereGeometry((x, y, z), r, periodic_boundary_condition=(1, 1, 1))
    tmp_carved = carve_geometry(amorph, geometry, side="in")
    print(f"tmp carved: {tmp_carved}")


amorph.write("amorph.data", format="lammps-data")

# Pack water into voids
water = pack_water(nummol=10000, atoms=amorph)
system = water + amorph
water.write("water.data", format="lammps-data")
system.write("system.data", format="lammps-data")
