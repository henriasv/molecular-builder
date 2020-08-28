"""Here, we demonstrate how to carve out a procedural surface
"""

from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import ProceduralSurfaceGeometry

# Create bulk of beta-cristobalite
atoms = create_bulk_crystal("beta_cristobalite", [120, 120, 120])

# Carve out geometry from beta-cristobalite
geometry = ProceduralSurfaceGeometry(center=(60,60,5),
                                     normal=(1,0,0),
                                     length=(120,120,10),
                                     method='simplex')
num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_with_procedural_surface.data", format="lammps-data")
carved.write(f"carved_surface.data", format="lammps-data")
