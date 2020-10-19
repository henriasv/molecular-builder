from molecular_builder.core import create_bulk_crystal, write, carve_geometry
from molecular_builder.geometry import NotchGeometry
import ase.cell 

L = [100, 100, 100]

quartz_structure = create_bulk_crystal("alpha_quartz", L) 
cell = quartz_structure.get_cell(); cell[0,1] = 0; cell[1,0] = 0; cell[0,2] = 0; cell[2,0] = 0; cell[1,2] = 0; cell[2,1] = 0
quartz_structure.set_cell(cell); quartz_structure.wrap()

crack = NotchGeometry([0,L[1]/2,L[1]/2], [20,0,0], [0,0,10])

carve_geometry(quartz_structure, crack, side="in")
write(quartz_structure, "alpha_quartz_with_crack.png",  camera_dir=[0, 1, 0], viewport_type="perspective")