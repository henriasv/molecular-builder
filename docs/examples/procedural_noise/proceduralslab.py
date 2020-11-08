from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import ProceduralSlabGeometry

atoms = create_bulk_crystal("beta_cristobalite", [100, 100, 50])

# Geometry-object for carving out a structure using Simplex/Perlin noise
slab_geometry = ProceduralSlabGeometry(point=(50, 50, 25),
                         normal=(0, 0, 1),
                         thickness=50,
                         octaves=1,
                         method='simplex',
                         )

num_carved, carved = carve_geometry(atoms, slab_geometry, side="in", return_carved=True)

atoms.write("block_with_procedural_slab.data", format='lammps-data')
carved.write("carved_slab.data", format='lammps-data')
