from molecular_builder import create_bulk_crystal, write
atoms = create_bulk_crystal("beta_cristobalite", size=[20,20,20])
write(atoms, "beta_cristobalite.data")
write(atoms, "beta_cristobalite.png", viewport_type="orthogonal")