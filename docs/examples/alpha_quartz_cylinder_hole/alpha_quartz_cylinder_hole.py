from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import CylinderGeometry

atoms = create_bulk_crystal("alpha_quartz", [60, 60, 60])

geometry = CylinderGeometry([30,30,30], 20, 100, orientation=[1,0,0])

carve_geometry(atoms, geometry, side="in")

write(atoms, "alpha_quartz_cylinder_hole.data")
write(atoms, "alpha_quartz_cylinder_hole.png", camera_dir=[3, 1, -1])