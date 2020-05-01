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
    def distance_point_line(n, q, p):
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
    def distance_point_plane(n, q, p):
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
        n = np.atleast_2d(n)    # Ensure n is 2d
        return np.abs(np.einsum('ik,jk->ij', p - q, n))

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
    orientation : nested list / ndarray_like
        orientation of cylinder, given as a vector pointing along the cylinder
        Random orientation by default
    kwargs : 
        properties
    """
    
    def __init__(self, center, length, orientation=[], **kwargs):
        super().__init__(**kwargs)
        assert len(center) == len(length), \
                 ("center and length need to have equal shapes")
        self.center = np.array(center)
        self.length = np.array(length) / 2
        
        # Set coordinate according to orientation
        if len(orientation) == 0:
            #orientation.append(np.random.randn(len(center)))
            orientation = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        if len(orientation) == 1:
            n_x = np.array(orientation[0])
            n_y = np.random.randn(len(center))
            n_y -= n_y.dot(n_x) * n_x
            orientation.append(n_y)
        if len(orientation) == 2:
            orientation.append(np.cross(orientation[0], orientation[1]))
        orientation = np.array(orientation, dtype=float) 
        self.orientation = orientation / np.linalg.norm(orientation, axis=1)
         
    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc) 
        indices = np.all((self.distance_point_plane(self.orientation, self.center, positions) <= self.length), axis=1)
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
            orientation = np.array(orientation, dtype=float)
            self.orientation = orientation / np.linalg.norm(orientation)
            
    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)
        
        indices = (self.distance_point_line(self.orientation, self.center, positions) <= self.radius) & \
                  (self.distance_point_plane(self.orientation, self.center, positions).flatten() <= self.length)
        return indices

class BerkovichGeometry(Geometry):
    def __init__(self, tip, axis=[0,0,-1], angle=np.radians(65.27)):
        self.indenter_angle = angle
        self.tip = np.asarray(tip)
        self.axis = np.asarray(axis)
        self.plane_directions = []
        self._create_plane_directions()

    def _create_plane_directions(self):
        xy_angles = [0, np.radians(120), np.radians(240)]
        for xy_angle in xy_angles:
            z_component = np.cos(np.pi/2-self.indenter_angle)
            xy_component = np.sin(np.pi/2-self.indenter_angle)
            self.plane_directions.append(np.asarray([
                                          xy_component*np.cos(xy_angle),
                                          xy_component*np.sin(xy_angle),
                                          z_component
                                          ]))

    def __call__(self, atoms):
        positions = atoms.get_positions()
        rel_pos = positions-self.tip
        is_inside_candidate1 = np.dot(rel_pos, self.plane_directions[0]) > 0
        is_inside_candidate2 = np.dot(rel_pos, self.plane_directions[1]) > 0
        is_inside_candidate3 = np.dot(rel_pos, self.plane_directions[2]) > 0
        is_inside = np.logical_and(np.logical_and(is_inside_candidate1, is_inside_candidate2), is_inside_candidate3)
        return is_inside
