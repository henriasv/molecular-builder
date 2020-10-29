from molecular_builder import create_bulk_crystal, write
atoms = create_bulk_crystal("alpha_quartz", size=[20,20,20])
write(atoms, "alpha_quartz.data")
write(atoms, "alpha_quartz.png", viewport_type="orthogonal")