from molecular_builder import write, create_bulk_ice

atoms = create_bulk_ice("Ih", [3, 3, 3])

write(atoms, "ice.data")
write(atoms, "ice.png") 