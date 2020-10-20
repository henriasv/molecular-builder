Brucite block in water
==========================

We make some includes 

.. literalinclude:: brucite_in_water.py
    :lines: 1-3 

Create a (too large) block of brucite_in_water

.. literalinclude:: brucite_in_water.py
    :lines: 5-6 

Since `create_bulk_crystal` makes a block of only whole unit cells, the lengths are probably not exactly `L`, therefore we update `L`

.. literalinclude:: brucite_in_water.py
    :lines: 8 

We then create a geometry and carve out only part of the brucite block. 

.. literalinclude:: brucite_in_water.py
    :lines: 10-14

Compute volumes and pack the right amount of water around the brucite block. 

.. literalinclude:: brucite_in_water.py
    :lines: 16-18

Finally, output the system as a LAMMPS data file. 

.. literalinclude:: brucite_in_water.py
    :lines: 20

In addition, we create two images of the system, one with a perspective projection and one with an orthogonal projection.

.. literalinclude:: brucite_in_water.py
    :lines: 21-22

.. figure:: brucite_in_water_perspective.png
    
    Brucite immersed in water. Perspective projection.

The same image in orthogonal 

.. figure:: brucite_in_water_orthogonal.png
    
    Brucite immersed in water. Orthogonal projection.

The complete script: 

.. literalinclude:: brucite_in_water.py