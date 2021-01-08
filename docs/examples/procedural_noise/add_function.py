from ase.io import read
from molecular_builder import carve_geometry, write
from molecular_builder.geometry import ProceduralSurfaceGeometry

atoms = read("block.data", format="lammps-data", style="molecular")


# Define function to add
def f(x, y):
    if x < 100:
        return - x/100
    else:
        return -1 + (x-100)/100


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
                                     f=f
                                     )
num_carved = carve_geometry(atoms, geometry, side="in")

write(atoms, "add_function.data")
write(atoms, "add_function.png", camera_dir=[0, 1, -1])
