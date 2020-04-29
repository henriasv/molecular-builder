from molecular_builder import create_bulk_crystal, carve_geometry, fetch_prepared_system
from molecular_builder.geometry import BerkovichGeometry, CylinderGeometry

slab = fetch_prepared_system("vashishta_1990_like_amorphous_silica/quench_950K")
slab.cell[2,2] = 80 # Expand cell in z direction to fit indenter 

indenter = create_bulk_crystal("diamond", (144, 144, 80), round="down")
carve_geometry(indenter, BerkovichGeometry((75, 75, 40)), side="out")
carve_geometry(indenter, CylinderGeometry((75, 75, 50), 60, 200, orientation=(0,0,1)), side="out")

atoms = slab+indenter

atoms.write("indenter_and_glass.data", format="lammps-data")
