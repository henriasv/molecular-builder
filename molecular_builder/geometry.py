import numpy as np
from ase import Atom

class Geometry:
    """Base class for geometries."""
    def __init__(self, periodic_boundary_condition = (False, False, False), minimum_image_convention=True):
        self.minimum_image_convention = minimum_image_convention
        self.periodic_boundary_condition = periodic_boundary_condition
        pass

    def __call__(self, atoms):
        """The empty geometry. False because we define no particle to be in the dummy geometry"""
        return np.zeros(len(atoms), dtype=np.bool)

    @staticmethod
    def distance_point_line(vec, point_line, point_ext):
        """Returns the (shortest) distance between a line parallel to
        a normal vector 'vec' through point 'point_line' and an external
        point 'point_ext'.

        :param vec: unit vector parallel to line
        :type vec: ndarray
        :param point_line: point on line
        :type point_line: ndarray
        :param point_ext: external points
        :type point_ext: ndarray

        """
        return np.linalg.norm(np.cross(vec, point_ext - point_line), axis=1)

    @staticmethod
    def distance_point_plane(vec, point_plane, point_ext):
        """Returns the (shortest) distance between a plane with normal vector
        'vec' through point 'point_plane' and a point 'point_ext'.

        :param vec: normal vector of plane
        :type vec: ndarray
        :param point_plane: point on line
        :type point_plane: ndarray
        :param point_ext: external point(s)
        :type point_ext: ndarray
        """
        vec = np.atleast_2d(vec)    # Ensure n is 2d
        return np.abs(np.einsum('ik,jk->ij', point_ext - point_plane, vec))

    @staticmethod
    def vec_and_point_to_plane(vec, point):
        """Returns the (unique) plane, given a normal vector 'vec' and a
        point 'point' in the plane. ax + by + cz - d = 0

        :param vec: normal vector of plane
        :type vec: ndarray
        :param point: point in plane
        :type point: ndarray
        """
        return np.array((*vec, np.dot(vec,point)))

    @staticmethod
    def cell2planes(cell, pbc):
        # 3 planes intersect the origin by ase design.
        a = cell[0]
        b = cell[1]
        c = cell[2]

        n1 = np.cross(a, b)
        n2 = np.cross(c, a)
        n3 = np.cross(b, c)

        #n1 = n1/np.dot(n1, n1)
        #n2 = n2/np.dot(n2, n2)
        #n3 = n3/np.dot(n3, n3)

        origin = np.array([0,0,0])+pbc/2
        top = (a+b+c)-pbc/2

        plane1 = Geometry.vec_and_point_to_plane(n1, origin)
        plane2 = Geometry.vec_and_point_to_plane(n2, origin)
        plane3 = Geometry.vec_and_point_to_plane(n3, origin)
        plane4 = Geometry.vec_and_point_to_plane(-n1, top)
        plane5 = Geometry.vec_and_point_to_plane(-n2, top)
        plane6 = Geometry.vec_and_point_to_plane(-n3, top)

        return [plane1, plane2, plane3, plane4, plane5, plane6]

    def packmol_structure(self, number, side):
        """ Make structure.

        :param number: Number of water molecules
        :type number: int
        :param side: Pack water inside/outside of geometry
        :type side: str
        :returns: String with information about the structure
        """
        structure = f"structure water.pdb\n"
        structure += f"  number {number}\n"
        structure += f"  {side} {self.__repr__()} "
        for param in self.params:
            structure += f"{param} "
        structure += "\nend structure\n"
        return structure

class PlaneBoundTriclinicGeometry(Geometry):
    def __init__(self, cell, pbc=0.0):
        self.planes = self.cell2planes(cell, pbc)
        self.ll_corner = [0,0,0]
        a = cell[0,:]
        b = cell[1,:]
        c = cell[2,:]
        self.ur_corner = a+b+c

    def packmol_structure(self, number, side):
        """ Make structure.

        :param number: Number of water molecules
        :type number: int
        :param side: Pack water inside/outside of geometry
        :type side: str
        :returns: String with information about the structure
        """
        if side == "inside":
            side = "over"
        elif side == "outside":
            side = "below"
        structure = f"structure water.pdb\n"
        structure += f"  number {number}\n"
        for plane in self.planes:
            structure += f"  {side} plane "
            for param in plane:
                structure += f"{param} "
            structure += "\n"
        structure += "end structure\n"
        return structure

    def __call__(self, position):
        raise NotImplementedError

class SphereGeometry(Geometry):
    """Spherical geometry.

    :param center: Center of sphere
    :type center: array_like
    :param radius: radius of sphere
    :type length: float
    """
    def __init__(self, center, radius, **kwargs):
        super().__init__(**kwargs)
        self.center = center
        self.radius = radius
        self.radius_squared = radius**2
        self.params = list(self.center) + [radius]
        self.ll_corner = np.array(center) - radius
        self.ur_corner = np.array(center) + radius

    def __repr__(self):
        return 'sphere'

    def __call__(self, atoms):
        atoms.append(Atom(position=self.center))
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        distances = atoms.get_distances(-1, list(range(len(atoms)-1)), mic=self.minimum_image_convention)
        atoms.pop()
        atoms.set_pbc(tmp_pbc)
        indices = distances**2 < self.radius_squared
        return indices

class CubeGeometry(Geometry):
    """Cubic geometry.

    :param center: Center of cube
    :type center: array_like
    :param length: length of each side
    :type length: float
    """
    def __init__(self, center, length, **kwargs):
        super().__init__(**kwargs)
        self.center = np.array(center)
        self.length_half = length / 2
        self.params = list(self.center - self.length_half) + [length]
        self.ll_corner = self.center - self.length_half
        self.ur_corner = self.center + self.length_half

    def __repr__(self):
        return 'cube'

    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)
        indices = np.all((np.abs(self.distance_point_plane(np.eye(3), self.center, positions)) <= self.length_half), axis=1)
        return indices

class BoxGeometry(Geometry):
    """Box geometry.

    :param ll_corner: lower-left corner of box
    :type center: array_like
    :param length: Length of box in each direction
    :type length: array_like
    """
    def __init__(self, center, length, **kwargs):
        super().__init__(**kwargs)
        self.length = length
        self.length_half = np.array(length) / 2
        self.center = np.array(center)
        self.ll_corner = self.center-self.length_half
        self.ur_corner = self.center+self.length_half
        self.params = list(self.ll_corner) + list(self.ur_corner)

    def __repr__(self):
        return 'box'

    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)
        indices = np.all((np.abs(self.distance_point_plane(np.eye(3), self.center, positions)) <= self.length_half), axis=1)
        return indices

    def volume(self):
        return np.prod(self.length)

class BlockGeometry(Geometry):
    """This is a more flexible box geometry, where the angle

    :param center: the center point of the block
    :type center: array_like
    :param length: the spatial extent of the block in each direction.
    :type length: array_like
    :param orientation: orientation of block
    :type orientation: nested list / ndarray_like

    NB: Does not support pack_water and packmol
    NB: This geometry will be deprecated
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

    def __repr__(self):
        return 'block'

    def packmol_structure(self, number, side):
        """ Make structure.
        """
        raise NotImplementedError("BlockGeometry does not support pack_water")

    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)
        indices = np.all((np.abs(self.distance_point_plane(self.orientation, self.center, positions)) <= self.length), axis=1)
        return indices

class PlaneGeometry(Geometry):
    """Remove all particles on one side of one or more planes. Can be used to
    form any 3d polygon, among other geometries

    :param point: point on plane
    :type point: array_like
    :param normal: vector normal to plane
    :type normal: array_like
    """
    def __init__(self, point, normal, **kwargs):
        super().__init__(**kwargs)
        assert len(point) == len(normal), "Number of given points and normal vectors have to be equal"

        self.point = np.atleast_2d(point)
        normal = np.atleast_2d(normal)
        self.normal = normal / np.linalg.norm(normal, axis=1)[:, np.newaxis]

    def packmol_structure(self, number, side):
        """ Make structure.

        :param number: Number of water molecules
        :type number: int
        :param side: Pack water inside/outside of geometry
        :type side: str
        :returns: String with information about the structure
        """
        if side == "inside":
            side = "over"
        elif side == "outside":
            side = "below"

        ds = np.einsum('ij,ij->j', self.point, self.normal)

        structure += f"structure water.pdb\n"
        structure += f"  number {number}\n"
        for plane in range(len(self.normal)):
            a, b, c = self.normal[side]
            d = ds[side]
            structure += f"  {side} plane {a} {b} {c} {d} \n"
        structure += "end structure\n"
        return structure

    def __call__(self, atoms):
        positions = atoms.get_positions()
        indices = np.all(np.einsum('ijk,ik->ij', self.point[:, np.newaxis] - positions, self.normal) > 0, axis=0)
        return indices


class CylinderGeometry(Geometry):
    """Cylinder object.

    :param center: the center point of the cylinder
    :type center: array_like
    :param radius: cylinder radius
    :type radius: float
    :param length: cylinder length
    :type length: float
    :param orientation: orientation of cylinder, given as a vector pointing along the cylinder. Pointing in x-direction by default.
    :type orientation: array_like
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
        self.params = list(center) + list(self.orientation) + [radius, length]

    def __repr__(self):
        return 'cylinder'

    def __call__(self, atoms):
        tmp_pbc = atoms.get_pbc()
        atoms.set_pbc(self.periodic_boundary_condition)
        positions = atoms.get_positions()
        atoms.set_pbc(tmp_pbc)

        indices = (self.distance_point_line(self.orientation, self.center, positions) <= self.radius) & \
                  (self.distance_point_plane(self.orientation, self.center, positions).flatten() <= self.length)
        return indices

class BerkovichGeometry(Geometry):
    # TODO: Implement support for packmol through plane geometry
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

    def packmol_structure(self, number, side):
        """ Make structure.
        """
        raise NotImplementedError("BerkovichGeometry is not yet supported by pack_water")

    def __call__(self, atoms):
        positions = atoms.get_positions()
        rel_pos = positions-self.tip
        is_inside_candidate1 = np.dot(rel_pos, self.plane_directions[0]) > 0
        is_inside_candidate2 = np.dot(rel_pos, self.plane_directions[1]) > 0
        is_inside_candidate3 = np.dot(rel_pos, self.plane_directions[2]) > 0
        is_inside = np.logical_and(np.logical_and(is_inside_candidate1, is_inside_candidate2), is_inside_candidate3)
        return is_inside

class EllipsoidGeometry(Geometry):
    """ Ellipsoid geometry, satisfies the equation

    (x - x0)^2   (y - y0)^2   (z - z0)^2
    ---------- + ---------- + ---------- = d
        a^2          b^2          c^2

    :param center: center of ellipsoid (x0, y0, z0)
    :type center: array_like
    :param length_axes: length of each axis (a, b, c)
    :type length_axes: array_like
    :param d: scaling
    :type d: float
    """

    # TODO: Add orientation argument

    def __init__(self, center, length_axes, d, **kwargs):
        super().__init__(**kwargs)
        self.center = np.asarray(center)
        self.length = np.asarray(length_axes)
        self.length_sqrd = self.length**2
        self.d = d
        self.params = list(self.center) + list(self.length) + [self.d]
        self.ll_corner = self.center - self.length
        self.ur_corner = self.center + self.length

    def __repr__(self):
        return 'ellipsoid'

    def __call__(self, atoms):
        positions = atoms.get_positions()
        positions_shifted_sqrd = (positions - self.center)**2
        LHS = np.divide(positions_shifted_sqrd, self.length_sqrd).sum(axis=1)
        indices = (LHS <= self.d)
        return indices

class EllipticalCylinderGeometry(Geometry):
    """ Elliptical Cylinder

    :param center: center of elliptical cylinder
    :type center: array_like
    :param a: axes along x-axis
    :type a: float
    :param b: axes along y-axis
    :type b: float
    :param length: length of cylinder
    :type length: float
    :param orientation: which way the cylinder should point
    :type orientation: ndarray

    NB: This geometry is not supported by packmol or pack_water
    """

    # TODO: Fix orientation argument (two separate orientations)

    def __init__(self, center, a, b, length, orientation=None, **kwargs):
        super().__init__(**kwargs)
        self.center = np.asarray(center)
        self.a_sqrd, self.b_sqrd = a**2, b**2
        self.length_half = np.asarray(length) / 2

        if orientation is None:
            self.orientation = np.zeros_like(center)
            self.orientation[0] = 1
        else:
            orientation = np.array(orientation, dtype=float)
            self.orientation = orientation / np.linalg.norm(orientation)

    def packmol_structure(self, number, side):
        """ Make structure.
        """
        raise NotImplementedError("EllipticalCylinderGeometry is not supported by pack_water")

    def __call__(self, atoms):
        positions = atoms.get_positions()
        positions_shifted_sqrd = (positions - self.center)**2
        ellipse = positions_shifted_sqrd[:,0]/self.a_sqrd + positions_shifted_sqrd[:,1]/self.b_sqrd
        indices = (ellipse <= 1) & (self.distance_point_plane(self.orientation, self.center, positions).flatten() <= self.length_half)
        return indices



class ProceduralSurfaceGeometry(Geometry):
    """Creates procedural noise on a surface defined by a point, a normal
    vector and a thickness.
    """

    def __init__(self,
                 point,
                 normal,
                 thickness,
                 scale=100,
                 method='simplex',
                 f=lambda x, y, z: 0,
                 **kwargs):
        assert len(point) == len(normal), \
            "Number of given points and normal vectors have to be equal"

        from noise import snoise3, pnoise3
        if method == "simplex":
            self.noise = snoise3
        elif method == "perlin":
            self.noise = pnoise3

        self.point = np.atleast_2d(point)
        normal = np.atleast_2d(normal)
        self.normal = normal / np.linalg.norm(normal, axis=1)[:, np.newaxis]
        self.thickness_half = thickness / 2
        self.scale = scale
        self.f = f
        self.kwargs = kwargs

    def packmol_structure(self, number, side):
        """ Make structure.
        """
        raise NotImplementedError(
            "ProceduralNoiseSurface is not supported by pack_water")

    def __call__(self, atoms):
        positions = atoms.get_positions()
        # calculate distance from particles to plane defined by normal and center
        distances = self.distance_point_plane(self.normal, self.point, positions)
        # find the points on plane
        point_plane = positions + np.einsum('ij,kl->jkl', distances, self.normal)
        # a loop is actually faster than an all-numpy implementation
        # since pnoise3/snoise3 are written in C++
        noises = np.empty(distances.shape)
        for i in range(len(self.normal)):
            for j, point in enumerate(point_plane[i]):
                noises[j] = self.noise(*(point/self.scale), **self.kwargs) + \
                            self.f(*point)

        dist = np.einsum('ijk,ik->ij', self.point[:, np.newaxis] - positions, self.normal)
        noises = noises.flatten() * self.thickness_half
        indices = np.all(dist > noises, axis=0)
        return indices
