# Molecular builder 

This builder relies heavily on the atomic simulation environment (ase). I particular we use ase to generate crystals from space group descriptions. 

## Description
This package facilitates building molecular systems. It has functions for creating bulk crystals and for fetching periodic boxes of non-crystalline systems from the internet. After generating bulk crystals, the user may carve out specific geometries, and combine several elements to create a system. 

We also support adding water using packmol. This requires having packmol installed. 

## Installation 
```
pip install git+https://github.com/henriasv/molecular-builder 
```

## Usage
Simple use case for carving out a cylinder from a bulk alpha-quartz crystal. 
```python 
atoms = create_bulk_system("alpha_quartz", [100,100,100])
carved_atoms = carve_geometry(CylinderGeometry(...), side="in")
atoms.save("output.data")
```

More advanced use case for carving out randomly placed and saving the carved out atoms and the resulting system in different files. 
```python 
atoms = create_bulk_system("beta_quartz", [200,200,200])
num_spheres = 10
pos_scaling = 30
r_scaling = 3

geometries = [SphereGeometry(np.asarray([i,j,k])**pos_scaling, r*r_scaling) for i,j,k,r in np.random.randn(num_spheres,4)] 

carved_atoms = Atoms()
for geometry in geometries:
    carved_atoms.append(atoms.carve(geometry))

atoms.save("system.data")
carved_atoms.save("carved_atoms.data")
```

Use case for setting up a system with a diamond indenter above a slab of beta quartz. 
```python 
quartz_atoms = create_bulk_system("beta_quartz", [200,200,200])
diamond_atoms = create_bulk_system("diamond", [180, 180, 200])

carved_atoms_quartz = quarts_atoms.carve(PlaneGeometry([0,0,1], [0,0,20]))
carved_atoms_diamond = diamond.carve(BerkovitzGeometry([0,0,-1], [0,0,30]))

diamond_atoms.center_on(quartz_atoms, "xy")
atoms = quartz_atoms+diamond_atoms

atoms.save("indenter_and_slab.data")
```

Use case for getting a particular input bulk structure. The name of the structure contains an identifier, in this case `p3754` that can be traced back to the procedure for creating the structure. For amorphous silica, this will typically be the melting and annealing process. 
```python 
atoms = fetch_bulk_system("amorphous_silica_p3754", [100,100,100])
atoms.save("my_amorphous_silica.data")
```

## Decisions to be made 
- Support periodic boundary conditions? If we are to, we need to have an on/off switch for this in the carve function, since some geometries, like notches, will behave poorly with periodic boundary conditions. 
- When adding atoms abject, which cell takes precendent? How to combine the cells of two atoms objects? 

