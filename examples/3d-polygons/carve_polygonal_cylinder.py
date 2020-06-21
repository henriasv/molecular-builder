"""Here, we demonstrate how the PlaneGeometry can be used to carve polygons of
an arbitrary number of corners and a certain height, named polygonal cylinders.

The slope of the unit circle at point (sin(θ), cos(θ)) is (cos(θ), sin(θ)).
"""

from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import PlaneGeometry

# Declare polygon parameters
num_sides = 5
scaling = 30
center = (60, 60, 0)

# Create bulk of beta-cristobalite
atoms = create_bulk_crystal("beta_cristobalite", [120, 120, 120])

# Define the various planes (normal vector and point in plane)
import numpy as np
thetas = np.linspace(0, 2*np.pi, num_sides+1)[:-1]
normals = [np.sin(thetas), np.cos(thetas), np.zeros(num_sides)]
normals = np.asarray(normals).T
points = scaling * normals + center

# Carve out geometry from beta-cristobalite
geometry = PlaneGeometry(points, normals)
num_carved, carved = carve_geometry(atoms, geometry, side="in", return_carved=True)

atoms.write("block_with_carved_pentagon.data", format="lammps-data")
carved.write(f"carved_pentagon.data", format="lammps-data")
