from molecular_builder import create_bulk_crystal
atoms = create_bulk_crystal("beta_cristobalite", size=[20,20,20])
atoms.write("system.data", format="lammps-data")