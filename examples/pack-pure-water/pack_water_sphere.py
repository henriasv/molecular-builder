""" In this example, we pack 10,000 water molecules in a sphere of
radius 50 Å and center (50, 50, 50). 
"""

from molecular_builder import pack_water
from molecular_builder.geometry import SphereGeometry

water = pack_water(nummol=10000, geometry=SphereGeometry((50, 50, 50), 50))
water.write("water_sphere.data", format="lammps-data")
