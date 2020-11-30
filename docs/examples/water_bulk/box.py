from molecular_builder import pack_water, write
from molecular_builder.geometry import BoxGeometry

geometry = BoxGeometry((50, 25, 25), (100, 50, 50))
atoms = pack_water(nummol=2000, geometry=geometry)

write(atoms, "box.data")
write(atoms, "box.png", camera_dir=[1, 2, -1])
