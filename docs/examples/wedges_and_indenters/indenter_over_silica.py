from molecular_builder import create_bulk_crystal, carve_geometry, fetch_prepared_system, write
from molecular_builder.geometry import BerkovichGeometry, CylinderGeometry

slab = fetch_prepared_system("vashishta_1990_like_amorphous_silica/quench_950K", type_mapping=[(1, 14), (2, 8)])
print(slab.numbers)

slab.cell[2,2] = 80 # Expand cell in z direction to fit indenter 

indenter = create_bulk_crystal("diamond", (144, 144, 80), round="down")
carve_geometry(indenter, BerkovichGeometry((75, 75, 40)), side="out")
carve_geometry(indenter, CylinderGeometry((75, 75, 50), 60, 200, orientation=(0,0,1)), side="out")

atoms = slab+indenter

write(atoms,"indenter_over_silica.png", bond_specs=[("Si", "O", 1.9)], camera_dir=(2, 1.5, -0.4), atom_radii=[("Si", 0.2), ("O", 0.2)], size=(1280, 960))

