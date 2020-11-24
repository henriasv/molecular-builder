"""Here, we demonstrate how the ProceduralSurfaceGeometry class
can be used to carve out a two-level surface
"""
import numpy as np
from molecular_builder import create_bulk_crystal, carve_geometry
from molecular_builder.geometry import ProceduralSurfaceGeometry

# User-defined variables
xlen = 100      # desired x length of system
ylen = 100      # desired y length of system
zlen = 50       # desired z length of system

noise_thickness = 10
scale = 4
octaves = 2

# Define geometrical properties of alpha-quartz
a = 5.02778179  # length of triclinic unit cell in a-dir (and x-dir)
b = 5.02778179  # length of triclinic unit cell in b-dir
c = 5.51891759  # length of triclinic unit cell in c-dir (and z-dir)
angle = 60
h = b * np.sin(np.deg2rad(angle))  # length of triclinic unit cell in y-dir

numcellsx = -(-xlen // a)   # number of cells in x-dir, rounded up
numcellsy = -(-ylen // h)   # number of cells in y-dir, rounded up
numcellsz = -(-zlen // c)   # number of cells in z-dir, rounded up

lengthx = numcellsx * a     # actual length in x-dir
lengthy = numcellsy * h     # actual length in y-dir
lengthz = numcellsz * c     # actual length in z-dir

lengthnoise = np.asarray((lengthx, lengthy, lengthz))
scalenoise = lengthnoise / scale

# Create bulk of alpha-quartz
quartz = create_bulk_crystal("alpha_quartz", (xlen, ylen, zlen))

# Carve out geometry from alpha-quartz
surface_geometry = ProceduralSurfaceGeometry(point=(xlen/2, ylen/2, zlen/2),
                                             normal=(0, 0, 1),
                                             thickness=noise_thickness,
                                             octaves=octaves,
                                             scale=scalenoise,
                                             pbc=lengthnoise,
                                             angle=angle,
                                             method='perlin',
                                             threshold=0)
carve_geometry(quartz, surface_geometry, side="in")
quartz.write("quartz.data", format="lammps-data")
