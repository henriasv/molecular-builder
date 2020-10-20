Getting started
=============================
This section gives a quick introduction to some of the things you can do with molecular-builder. 

Creating bulk crystal structures
----------------------------------
The perhaps simples use case is to create a bulk structure of a given size. In this way, you leverage the database of crystsal structure specifications provided by molecular-builder. Creating a block of :math:`\beta`-cristobalite and saving it as a lammps data file may be done as follows:

.. literalinclude:: python/bulk_beta_cristobalite.py

The atoms variable now contains an `ase.Atoms` object, and the writing is the write method of that object. 

.. image:: python/beta_cristobalite.png


Similarly, one may create a block of :math:`\alpha`-quartz: 

.. literalinclude:: python/bulk_alpha_quartz.py

Since the unit cell of :math:`\alpha`-quartz is triclinic, the resulting structure may be triclinic, depending on whether the structure can be represented in an orhogonal cell. 

.. image:: python/alpha_quartz.png

The functionality for creating bulk structures uses the `ase.spacegroup` module, and simplifies the process by providing the parameters for the spacegroup command simply by providing a string with the name of the structure. A complete list of available structures can be found here @TODO. 

Carving out a sphere 
---------------------------------------


Carving out a block 
---------------------------------------


Carving out multiple structures 
---------------------------------------

