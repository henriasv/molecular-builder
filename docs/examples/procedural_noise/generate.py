from molecular_builder import create_bulk_crystal, carve_geometry, write
from molecular_builder.geometry import ProceduralSurfaceGeometry

# Create bulk of beta-cristobalite
atoms = create_bulk_crystal("beta_cristobalite", [200, 200, 50])

write(atoms, "block.data")
write(atoms, "block.png")
