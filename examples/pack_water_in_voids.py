import numpy as np
from molecular_builder import fetch_prepared_system, carve_geometry, pack_water
from molecular_builder.geometry import SphereGeometry, CubeGeometry

# Amorphous silica of size (357, 143, 143)
amorph = fetch_prepared_system("amorphous_silica_1")

num_spheres = 20
for sphere in range(num_spheres):
    i, j, k, l = np.random.uniform(size=4)
    x, y, z, r = i*357, j*143, k*143, l*30
    geometry = SphereGeometry([x, y, z], r, periodic_boundary_condition=(True, True, True))
    tmp_carved = carve_geometry(amorph, geometry, side="in")
    print(f"tmp carved: {tmp_carved}")

water = pack_water(10000, atoms=amorph) #, geometry=CubeGeometry([60,60,60], 130))
amorph.write("amorph.data", format="lammps-data")
water.write("out.data", format="lammps-data")
