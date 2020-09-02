from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import BoxGeometry

atoms = create_bulk_crystal("silicon_carbide_3c", [100, 100, 100])

geometry = BoxGeometry([50,50,50], [70,70,70])

carve_geometry(atoms, geometry, side="out")

write(atoms, "box.data")
write(atoms, "box.png", camera_dir=[3, 2, -1])