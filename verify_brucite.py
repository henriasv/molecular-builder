#!/usr/bin/env python3
"""Verify the fixed brucite has only ONE Mg layer."""

import matplotlib
matplotlib.use('Agg')
from molecular_builder import create_bulk_crystal
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

# Create single-layer brucite using the fix
base = create_bulk_crystal('brucite', [3.5, 3.5, 3])  # Z=3 for single layer
slab = base.repeat((3, 3, 1))

# Unwrap H atoms
positions = slab.get_positions()
symbols = slab.get_chemical_symbols()
cell_z = slab.get_cell()[2, 2]

for i, (pos, sym) in enumerate(zip(positions, symbols)):
    if sym == 'H' and pos[2] < cell_z * 0.3:
        positions[i, 2] += cell_z
        
slab.set_positions(positions)

# Add vacuum
cell = slab.get_cell()
cell[2, 2] += 15.0
slab.set_cell(cell)
slab.center(axis=2)

# Check Mg layers
pos = slab.get_positions()
syms = slab.get_chemical_symbols()
mg_z = pos[[s == 'Mg' for s in syms], 2]

import numpy as np
unique_mg_z = np.unique(np.round(mg_z, 1))

print(f"Single-layer brucite verification:")
print(f"  Total atoms: {len(slab)}")
print(f"  Mg atoms: {len(mg_z)}")
print(f"  Unique Mg Z-levels: {len(unique_mg_z)}")
print(f"  Mg Z-positions: {unique_mg_z}")

if len(unique_mg_z) == 1:
    print("  ✓ SUCCESS: Only ONE Mg layer!")
else:
    print(f"  ✗ FAIL: {len(unique_mg_z)} Mg layers detected")

# Visualize
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
plot_atoms(slab, ax, rotation=('90x,0y,0z'))
ax.set_title(f'Fixed Brucite - Single Layer\n({len(unique_mg_z)} Mg layer, {len(slab)} atoms)')
plt.tight_layout()
plt.savefig('/Users/henriasv/repos/molecular-builder/brucite_single_layer.png', dpi=150)
print("\nSaved to brucite_single_layer.png")
