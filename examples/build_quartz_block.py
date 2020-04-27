from molecular_builder import create_bulk_crystal

atoms = create_bulk_crystal("alpha_quartz", [50,100,200])
atoms.write("alpha_quartz.data", format="lammps-data")