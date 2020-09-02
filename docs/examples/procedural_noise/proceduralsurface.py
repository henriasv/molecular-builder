from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import ProceduralSurfaceGeometry

# Create bulk of beta-cristobalite
atoms = create_bulk_crystal("beta_cristobalite", [200, 200, 50])
write(atoms, "block.png")

# Carve out geometry from beta-cristobalite
geometry = ProceduralSurfaceGeometry(point=(100, 100, 40),
                                     normal=(0, 0, 1),
                                     thickness=20,
                                     octaves=1,
                                     method='simplex',
                                     )
num_carved, carved = carve_geometry(atoms, geometry, side="out", return_carved=True)

write(atoms, "block_with_procedural_surface.png")