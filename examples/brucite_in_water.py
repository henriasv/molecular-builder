from molecular_builder import create_bulk_crystal, carve_geometry, write, pack_water
from molecular_builder.geometry import BlockGeometry
from ase.constraints import FixBondLengths

atoms = create_bulk_crystal("brucite", [120, 120, 120])

geometry = BlockGeometry(center=[60,60,59.5], 
                         length=[41,150,95], 
                         orientation=[[1,0,0],[0,1,0],[0,0,1]],
                         periodic_boundary_condition=(True, True, True))

num_carved = carve_geometry(atoms, geometry, side="out")

water = pack_water(20000, atoms=atoms, pbc=2.0)
atoms += water

write(atoms, "brucite_block.data", bond_specs=("O", "H", 1.02))
write(atoms, "brucite_block.png", bond_specs=("O", "H", 1.02))
