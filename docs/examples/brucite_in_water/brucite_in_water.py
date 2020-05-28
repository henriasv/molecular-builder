from molecular_builder import create_bulk_crystal, carve_geometry, write, pack_water
from molecular_builder.geometry import BlockGeometry
from ase.constraints import FixBondLengths
import numpy as np 

L = np.array([80, 60, 80])

atoms = create_bulk_crystal("brucite", L)

L = atoms.cell.lengths()
volume = atoms.cell.volume

geometry = BlockGeometry(center=L/2, 
                         length=[L[0], L[1]/2-3, 4*L[2]/5+2], 
                         orientation=[[1,0,0],[0,1,0],[0,0,1]],
                         periodic_boundary_condition=(True, True, True))

water_volume = volume-geometry.volume() 

carve_geometry(atoms, geometry, side="out")
volume=water_volume
pack_water(atoms=atoms, volume=water_volume, pbc=2.0, tolerance=1.8)

write(atoms, "brucite_in_water.data", bond_specs=("O", "H", 1.02))
write(atoms, "brucite_in_water_perspective.png", bond_specs=("O", "H", 1.02), camera_dir=[1, 0, 0], viewport_type="perspective", size=(480, 640))
write(atoms, "brucite_in_water_orthogonal.png", bond_specs=("O", "H", 1.02), camera_dir=[1, 0, 0], viewport_type="orthogonal", size=(480, 640))
