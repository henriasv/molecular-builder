#!/usr/bin/env python3
"""Debug what create_bulk_crystal actually creates for brucite."""

from molecular_builder import create_bulk_crystal
import numpy as np

# Test different sizes
for size_z in [3, 4, 5, 6, 10]:
    print(f"\n{'='*60}")
    print(f"Testing create_bulk_crystal('brucite', [3.5, 3.5, {size_z}])")
    print(f"{'='*60}")
    
    brucite = create_bulk_crystal('brucite', [3.5, 3.5, size_z])
    
    print(f"Total atoms: {len(brucite)}")
    print(f"Cell: {brucite.get_cell().tolist()}")
    print(f"Cell Z-height: {brucite.get_cell()[2, 2]:.3f} Å")
    
    # Analyze by element
    positions = brucite.get_positions()
    symbols = brucite.get_chemical_symbols()
    
    for element in ['Mg', 'O', 'H']:
        z_pos = positions[[s == element for s in symbols], 2]
        if len(z_pos) > 0:
            unique_z = np.unique(np.round(z_pos, 2))
            print(f"  {element}: {len(z_pos)} atoms, unique Z-levels: {len(unique_z)}")
            print(f"       Z range: {z_pos.min():.3f} to {z_pos.max():.3f} Å")
            if element == 'Mg':
                print(f"       Unique Z positions: {sorted(unique_z)[:5]}...")  # Show first 5

print(f"\n{'='*60}")
print("CONCLUSION: Check if Z-size affects number of Mg layers")
