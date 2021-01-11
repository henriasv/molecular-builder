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
                                     seed=16591,
                                     threshold=0,
                                     repeat=True,
                                     angle=45
                                     )
num_carved = carve_geometry(atoms, geometry, side="out")

write(atoms, "screwed.data")
write(atoms, "screwed.png", camera_dir=[0, 0, -1])
