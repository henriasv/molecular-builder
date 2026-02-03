""" In this example, we pack 10,000 water molecules in a sphere of
radius 50 Å and center (50, 50, 50). 
"""

from molecular_builder import pack_water
from molecular_builder.geometry import SphereGeometry

# Volume of sphere R=50 is ~523,000 A^3. @ 0.033 density -> ~17,000 waters.
# 10,000 is actually fine here! (Density ~0.6 g/cm3). 
# But let's leave it as is, or explicit about it.
water = pack_water(nummol=10000, geometry=SphereGeometry((50, 50, 50), 50))
water.write("water_sphere.data", format="lammps-data")
