# molecular-builder
Builder for molecular systems


Built on atomic simulation environment (ase). 


```
atoms = create_bulk_system("alpha_quartz", [100,100,100])
carved_atoms = atoms.carve(CylinderGeometry(...))
atoms.save("output.data")
```

```
atoms = create_bulk_system("beta_quartz", [200,200,200])

```
