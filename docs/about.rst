About
==============
Molecular builder is a python package for creating input structures for molecular dynamics simulations. 

It is focused on creating bulk structures and carving out geometries. Therefore, `molecular-builder` comes with a database of crystal structures and a collection og geometries that can be used to create various structures in various shapes. 

A main goal of this package is to make the whole process of setting up an initial geometry for a  molecular dynamics simulation programmable, and thus reliable.

For example, creating an :math:`\alpha`-quartz block with a cylidrical hole in it can be done as follows: 

.. literalinclude:: examples/alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.py

This simple script produces a lammps data file and an image of the created system: 

.. image:: examples/alpha_quartz_cylinder_hole/alpha_quartz_cylinder_hole.png
