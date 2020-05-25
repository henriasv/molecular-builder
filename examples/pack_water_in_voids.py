import numpy as np
from molecular_builder import create_bulk_crystal, carve_geometry, pack_water
from molecular_builder.geometry import SphereGeometry, CubeGeometry

atoms = create_bulk_crystal("alpha_quartz", [120, 120, 120])

num_spheres = 10
for sphere in range(num_spheres):
    i, j, k, l = np.random.uniform(size=4)
    x, y, z, r = i*120, j*120, k*120, l*30
    geometry = SphereGeometry([x, y, z], r, periodic_boundary_condition=(True, True, True))
    tmp_carved = carve_geometry(atoms, geometry, side="in")
    print(f"tmp carved: {tmp_carved}")

water = pack_water(atoms, 1000, CubeGeometry([60,60,60], 130))
