from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import BoxGeometry

atoms = create_bulk_crystal("alpha_quartz", (120, 120, 120))

geometry = BoxGeometry(center=(40, 40, 40), length=(60, 60, 60))

num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_with_hole_block.data", format="lammps-data")
carved.write("carved_block.data", format="lammps-data")
