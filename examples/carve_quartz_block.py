from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import BlockGeometry

atoms = create_bulk_crystal("alpha_quartz", [120, 120, 120])

geometry = BlockGeometry([60,90,60], [80,80,80], periodic_boundary_condition=(True, True, True))

num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_with_hole_block.data", format="lammps-data")
carved.write("carved_block.data", format="lammps-data")
