import numpy as np 
from ase import Atom

class Geometry:
    def __init__(self, periodic_boundary_condition = (False, False, False), minimum_image_convention=True):
        self.minimum_image_convention = minimum_image_convention
        self.periodic_boundary_condition = periodic_boundary_condition
        pass
    
    def __call__(self, atoms):
        """The empty geometry. False because we define no particle to be in the dummy geometry"""
        return np.zeros(len(atoms), dtype=np.bool) 

class SphereGeometry(Geometry):
    def __init__(self, center, radius, **kwargs):
        super().__init__(**kwargs)
        self.center = center
        self.radius = radius 
        self.radius_squared = radius**2
         
    def __call__(self, atoms):
        atoms.append(Atom(position=self.center))
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        distances = atoms.get_distances(-1, list(range(len(atoms)-1)), mic=self.minimum_image_convention)
        print(len(distances))
        atoms.pop()
        atoms.set_pbc(tmp_pbc)
        indices = distances**2 < self.radius_squared
        print(indices)
        print(len(indices))
        return indices
        
class BlockGeometry(Geometry):
    """ Block geometry gives the indices of a given block.
    
    Parameters
    ----------
    center : array_like
        the center point of the block (size needs to be equal to number of 
        spatial dimensions)
    length : array_like
        the spatial extent of the block in each direction. 
    kwargs : 
        properties
    """
    
    def __init__(self, center, length, **kwargs):
        super().__init__(**kwargs)
        assert len(center) == len(length), \
                 ("center and length need to have equal shapes")
        self.center = np.array(center)
        self.length = np.array(length)
        self.llcorner = self.center - self.length/2   # Lower left corner
        self.urcorner = self.center + self.length/2   # Upper right corner

         
    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)
        indices = np.all((self.llcorner <= positions) & 
                         (positions <= self.urcorner), axis=1)
        return indices

