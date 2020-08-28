"""Here, we demonstrate how to carve out a procedural surface
"""

from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import ProceduralSurfaceGeometry

# Create bulk of beta-cristobalite
atoms = create_bulk_crystal("beta_cristobalite", [50, 200, 200])
atoms.write("system.data", format="lammps-data")

# Carve out geometry from beta-cristobalite
geometry = ProceduralSurfaceGeometry(point=(40, 100, 100),
                                     normal=(1, 0, 0),
                                     thickness=20,
                                     octaves=1,
                                     method='simplex')
num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_with_procedural_surface.data", format="lammps-data")
carved.write(f"carved_surface.data", format="lammps-data")
