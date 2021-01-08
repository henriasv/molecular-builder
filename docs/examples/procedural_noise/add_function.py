from math import exp
from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import ProceduralSurfaceGeometry

# Create bulk of beta-cristobalite
atoms = create_bulk_crystal("beta_cristobalite", [200, 200, 50])
write(atoms, "block.png")


# Define function to add
def f(x, y):
    val = exp(-((x-100)**2 + (y-100)**2)/1000)
    print(x, y, val)
    return val


# Carve out geometry from beta-cristobalite
geometry = ProceduralSurfaceGeometry(point=(100, 100, 30),
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
write(atoms, "add_function.png", camera_dir=[0, 0, -1])
