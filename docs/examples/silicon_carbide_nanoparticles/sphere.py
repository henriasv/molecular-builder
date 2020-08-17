from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import SphereGeometry

atoms = create_bulk_crystal("silicon_carbide_3c", [100, 100, 100])

geometry = SphereGeometry((50,50,50), 40)

carve_geometry(atoms, geometry, side="out")

write(atoms, "sphere.data")
write(atoms, "sphere.png", camera_dir=[3, 1, -1])