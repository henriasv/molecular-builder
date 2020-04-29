from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import CylinderGeometry

atoms = create_bulk_crystal("alpha_quartz", [120, 120, 120])

geometry = CylinderGeometry([60,60,80], 20, 90, orientation=[0.5,1,1])

num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_with_hole_cylinder.data", format="lammps-data")
carved.write("carved_cylinder.data", format="lammps-data")
