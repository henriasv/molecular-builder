Water bulk
=================
Bulks of pure water can be packed using the `pack_water` method implemented in molecular-builder. The functionality is based on `Packmol <http://m3g.iqm.unicamp.br/packmol/home.shtml>`_. Make sure that Packmol installed before using the function.

First, a geometry where the water should be packed has to be defined. Most of the geometries available in molecular-builder (but not all of them) are supported by Packmol, including:

- `SphereGeometry`
- `CubeGeometry`
- `BoxGeometry`
- `PlaneGeometry`
- `CylinderGeometry`
- `EllipsoidGeometry`
- `OctahedronGeometry`
- `DodecahedronGeometry`

Then, the geometry is sent into `pack_water` using the `geometry` argument. Below, an example where 2000 water molecules are packed into a box of dimensions (100Å, 50Å, 50Å) is given.

.. literalinclude:: box.py

.. figure:: box.png
