![build docs](https://github.com/henriasv/molecular-builder/workflows/build%20docs/badge.svg)

## Molecular builder
This package facilitates building molecular systems. It has functions for creating bulk crystals and for fetching periodic boxes of non-crystalline systems from the internet. After generating bulk crystals, the user may carve out specific geometries, and combine several elements to create a system. 

We also support adding water using packmol. This requires having packmol installed. 

# Dependencies 
This molecular builder relies heavily on the atomic simulation environment (ase) and ovito. Ase is used to generate crystals from space group descriptions, ovito is used to create bonds and to output molecular systems in the lammps data format. 

This molecular builder is made with LAMMPS in mind. Usually, saving to input files for other simulators will work just fine, but anywhere a choice has to be made, convenience when preparing systems for LAMMPS will be prioritized. 


## Installation 
```
pip install molecular-builder 
```

## Usage
Simple use case for creating a bulk alpha-quartz crystal. 
```python 
from molecular_builder import create_bulk_crystal
atoms = create_bulk_crystal("alpha_quartz", [50,100,200])
write(atoms, "alpha_quartz.data)
```

Please refer to the [documentation](https://henriasv.github.io/molecular-builder) for more usage examples. 