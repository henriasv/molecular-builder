from molecular_builder import create_bulk_crystal

atoms = create_bulk_crystal("silicon_carbide_3c_110", (1, 1, 1))
atoms.write("sic_110.data", format="lammps-data")
