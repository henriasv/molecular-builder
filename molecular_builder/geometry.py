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
        
    @staticmethod
    def distancePointLine(n, q, p):
        """ Returns the (shortest) distance between a line parallel to
        a normal vector n through point q and a point p.
        
        Parameters
        ----------
        n : ndarray
            unit vector parallel to line
        q : ndarray
            point on line
        p : ndarray
            external points
        """
        return np.linalg.norm(np.cross(n, p - q), axis=1)
        
    @staticmethod
    def distancePointPlane(n, q, p):
        """ Returns the (shortest) distance between a plane with normal vector 
        n through point q and a point p.
        
        Parameters
        ----------
        n : ndarray
            unit vector normal to plane
        q : ndarray
            point in plane
        p : ndarray
            external points
        """
        return np.abs((p - q).dot(n))

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
        atoms.pop()
        atoms.set_pbc(tmp_pbc)
        indices = distances**2 < self.radius_squared
        return indices
        
class BlockGeometry(Geometry):
    """ Block object.
    
    Parameters
    ----------
    center : array_like
        the center point of the block
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
        
class CylinderGeometry(Geometry):
    """ Cylinder object.
    
    Parameters
    ----------
    center : array_like
        the center point of the cylinder
    radius : float
        cylinder radius
    length : float
        cylinder length
    orientation : array_like
        orientation of cylinder, given as a vector pointing along the cylinder
        Pointing in x-direction by default.
    """
    
    def __init__(self, center, radius, length, orientation=None, **kwargs):
        super().__init__(**kwargs)
        self.center = np.array(center)
        self.radius = radius
        self.length = length / 2
        if orientation is None:
            self.orientation = np.zeros_like(center)
            self.orientation[0] = 1
        else:
            self.orientation = np.array(orientation) / np.linalg.norm(np.array(orientation))
            
            
    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)
        indices = (self.distancePointLine(self.orientation, self.center, positions) <= self.radius) & \
                  (self.distancePointPlane(self.orientation, self.center, positions) <= self.length)
        return indices

