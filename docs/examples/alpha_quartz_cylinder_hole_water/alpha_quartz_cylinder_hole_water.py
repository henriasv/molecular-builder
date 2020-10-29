from molecular_builder import create_bulk_crystal, carve_geometry, write, pack_water
from molecular_builder.geometry import CylinderGeometry
import numpy as np

atoms = create_bulk_crystal("alpha_quartz", [60, 60, 60])

r = 25
l = atoms.cell[0][0]

geometry = CylinderGeometry([30,30,30], r, 100, orientation=[1,0,0])

water_volume = np.pi*r**2*l

carve_geometry(atoms, geometry, side="in")

pack_water(atoms, volume=water_volume, pbc=4.0, tolerance=1.5, density=0.8)

write(atoms, "alpha_quartz_cylinder_hole_water.data")
write(atoms, "alpha_quartz_cylinder_hole_water.png", camera_dir=[3, 1, -1])