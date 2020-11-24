from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import OctahedronGeometry

atoms = create_bulk_crystal("silicon_carbide_3c", [100, 100, 100])

geometry = OctahedronGeometry(30, (50, 50, 50))

carve_geometry(atoms, geometry, side="out")

write(atoms, "octahedron.data")
write(atoms, "octahedron.png", camera_dir=[3, 1, -1])
