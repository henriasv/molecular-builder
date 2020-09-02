from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import CylinderGeometry

atoms = create_bulk_crystal("silicon_carbide_3c", [100, 100, 100])

geometry = CylinderGeometry([50,50,50], 20, 100, orientation=[1,0,1])

carve_geometry(atoms, geometry, side="out")

write(atoms, "cylinder.data")
write(atoms, "cylinder.png", camera_dir=[3, 2, -1])