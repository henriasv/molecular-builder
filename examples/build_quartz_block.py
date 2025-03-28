from molecular_builder import create_bulk_crystal, write

atoms = create_bulk_crystal("alpha_quartz", [50,100,50])

write(atoms, "alpha_quartz_mb.data", bond_specs=[("Si", "O", 2.8)], specorder=["O", "Si"])
write(atoms, "alpha_quartz_mb.png", bond_specs=[("Si", "O", 2.8)])