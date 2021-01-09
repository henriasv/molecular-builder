from ase.io import read
from molecular_builder import carve_geometry, write
from molecular_builder.geometry import ProceduralSurfaceGeometry

atoms = read("block.data", format="lammps-data", style="molecular")

# Carve out geometry from beta-cristobalite
geometry = ProceduralSurfaceGeometry(point=(100, 100, 40),
                                     normal=(0, 0, 1),
                                     thickness=20,
                                     method='simplex',
                                     scale=3,
                                     octaves=2,
                                     seed=16591
                                     )
num_carved = carve_geometry(atoms, geometry, side="out")

write(atoms, "parameters.data")
write(atoms, "parameters.png", camera_dir=[2, 1, -1])
