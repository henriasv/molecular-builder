#!/usr/bin/env python3
"""Test the fixed brucite crystal definition."""

from molecular_builder import create_bulk_crystal
import numpy as np

print("Testing fixed brucite crystal definition...")
print("="*60)

# Create single-layer brucite
brucite = create_bulk_crystal('brucite', [3.5, 3.5, 3])

pos = brucite.get_positions()
syms = brucite.get_chemical_symbols()
cell_z = brucite.get_cell()[2, 2]

print(f"\nCell Z-height: {cell_z:.3f} Å")
print(f"Total atoms: {len(brucite)}")

# Analyze positions
for element in ['Mg', 'O', 'H']:
    z = pos[[s == element for s in syms], 2]
    if len(z) > 0:
        print(f"\n{element}:")
        print(f"  Count: {len(z)}")
        print(f"  Z range: {z.min():.3f} to {z.max():.3f} Å")
        print(f"  Z center: {(z.min() + z.max())/2:.3f} Å")

# Check if Mg is centered
mg_z = pos[[s == 'Mg' for s in syms], 2]
mg_center = (mg_z.min() + mg_z.max()) / 2
cell_center = cell_z / 2

print(f"\n{'='*60}")
print(f"Centering check:")
print(f"  Mg layer center: {mg_center:.3f} Å")
print(f"  Cell center: {cell_center:.3f} Å")
print(f"  Offset: {abs(mg_center - cell_center):.3f} Å")

if abs(mg_center - cell_center) < 0.1:
    print("  ✓ Mg is CENTERED in cell!")
else:
    print("  ✗ Mg is NOT centered")

# Check OH symmetry
o_z = pos[[s == 'O' for s in syms], 2]
h_z = pos[[s == 'H' for s in syms], 2]

print(f"\nOH group symmetry:")
print(f"  O below Mg: {np.sum(o_z < mg_center)} atoms")
print(f"  O above Mg: {np.sum(o_z > mg_center)} atoms")
print(f"  H below Mg: {np.sum(h_z < mg_center)} atoms")
print(f"  H above Mg: {np.sum(h_z > mg_center)} atoms")
