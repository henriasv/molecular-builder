
Cracks 
==================

To create a v-shaped defect in a material, the Notch geometry can be used. In this example we create a slab of :math:`\alpha`-quartz and carve out a crack made with the Notch geometry.

First, we import the necessary packages: 

.. literalinclude:: crack_in_quartz.py
    :lines: 1-3

Then, we create our slab of :math:`\alpha`-quartz. We use block dimensions of (100Å, 50Å, 100Å)

.. literalinclude:: crack_in_quartz.py
    :lines: 5-9

The code below will create and carve out a crack in the :math:`\alpha`-quartz using the Notch geometry: 

.. literalinclude:: crack_in_quartz.py
    :lines: 11-12

The resulting system looks like this:

.. figure:: alpha_quartz_with_crack.png
    
    An crack in a slab of :math:`\alpha`-quartz.

Full code of this example:

.. literalinclude:: crack_in_quartz.py
