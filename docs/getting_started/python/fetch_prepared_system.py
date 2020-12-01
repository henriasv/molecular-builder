from molecular_builder import fetch_prepared_system, write

atoms = fetch_prepared_system("vashishta_1990_like_amorphous_silica_quench_300K", type_mapping=[[1, 14], [2, 8]])
write(atoms, "amorphous_silica.png")