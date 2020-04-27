from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import SphereGeometry

atoms = create_bulk_crystal("alpha_quartz", [120, 120, 120])

geometry = SphereGeometry([60,90,60], 40, periodic_boundary_condition=(True, True, True))

block_with_hole = carve_geometry(atoms, geometry, side="in")
carved = carve_geometry(atoms, geometry, side="out")

block_with_hole.write("block_with_hole.data", format="lammps-data")
carved.write("carved.data", format="lammps-data")
