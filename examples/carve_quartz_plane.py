from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import PlaneGeometry

atoms = create_bulk_crystal("alpha_quartz", [120, 120, 120])

geometry = PlaneGeometry([60,60,60], [1,0,0], periodic_boundary_condition=(True, True, True))

num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_inside_plane.data", format="lammps-data")
carved.write("block_outside_plane.data", format="lammps-data")
