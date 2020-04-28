from molecular_builder import create_bulk_crystal, carve_geometry, fetch_prepared_system
from molecular_builder.geometry import SphereGeometry
import numpy as np 

atoms = fetch_prepared_system("amorphous_silica_1")

num_spheres = 20

geometries = [SphereGeometry(np.asarray([i*357,j*143,k*143]), 
                                r*30,
                                periodic_boundary_condition=(True, True, True)) for i,j,k,r in np.random.uniform(size=(num_spheres,4))] 

num_carved = 0
for geometry in geometries:
    tmp_carved = carve_geometry(atoms, geometry, side="in")
    print(f"tmp carved: {tmp_carved}")
    num_carved += tmp_carved

print(f"Carved out {num_carved} atoms")

atoms.write("block_with_holes.data", format="lammps-data")
