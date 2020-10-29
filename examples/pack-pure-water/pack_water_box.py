""" In this example, we pack 10,000 water molecules in a box with center
(50, 25, 25) and length (100Å, 50Å, 50Å). 
"""

from molecular_builder import pack_water
from molecular_builder.geometry import BoxGeometry

water = pack_water(nummol=10000, geometry=BoxGeometry((50, 25, 25), (100, 50, 50)))
water.write("water_box.data", format="lammps-data")
