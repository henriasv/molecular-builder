#!/usr/bin/env python3
"""
Benchmark for molecular-builder native water packing.

Creates periclase (MgO) and brucite (Mg(OH)2) slabs at different sizes,
packs water using the native method, and saves all results to a single extxyz file.
"""

import time
import numpy as np
from ase import Atoms
from ase.io import write
from molecular_builder import create_bulk_crystal, pack_water


def create_periclase_slab(size):
    """Create MgO slab using ASE bulk."""
    from ase.build import bulk
    print(f"  Creating periclase slab (~{size}A)...")
    
    # Create bulk MgO (rocksalt structure)
    n = max(2, int(round(size / 4.21)))
    mgo = bulk('MgO', 'rocksalt', a=4.21, cubic=True).repeat((n, n, max(1, int(n*0.4))))
    
    # Add vacuum for water
    cell = mgo.get_cell()
    cell[2, 2] += max(10.0, size * 0.5)
    mgo.set_cell(cell)
    mgo.center(axis=2)
    
    return mgo


def create_brucite_slab(size):
    """Create single-layer Brucite slab (now properly centered in crystal definition)."""
    print(f"  Creating brucite slab (~{size}A)...")
    
    # Z=3 gives single Mg layer, properly centered with OH groups above and below
    base = create_bulk_crystal('brucite', [3.5, 3.5, 3])
    
    # Repeat only in X and Y to maintain single layer
    n_xy = max(2, int(round(size / 3.142)))
    slab = base.repeat((n_xy, n_xy, 1))
    
    # Add vacuum for water
    cell = slab.get_cell()
    cell[2, 2] += max(15.0, size * 0.6)
    slab.set_cell(cell)
    slab.center(axis=2)
    
    return slab


def run_benchmark():
    """Run benchmark on different sizes."""
    sizes = [10, 15, 20]
    systems = ['Periclase', 'Brucite']
    density = 1.0
    
    all_structures = []
    results = []
    
    print("="*60)
    print("Molecular-Builder Native Water Packing Benchmark")
    print("="*60)
    
    for system_name in systems:
        for size in sizes:
            print(f"\n{system_name} (size ~{size}A):")
            
            # Create system
            if system_name == 'Periclase':
                atoms = create_periclase_slab(size)
            else:
                atoms = create_brucite_slab(size)
            
            n_substrate = len(atoms)
            
            # Pack water with native method
            print(f"  Packing water (native method, density={density} g/cm3)...")
            start = time.time()
            
            water = pack_water(
                atoms,
                method="native",
                density=density,
                seed=42,
                pairwise_distances={('O', 'Mg'): 2.0}
            )
            
            elapsed = time.time() - start
            
            # Results
            n_waters = len(water) // 3
            total_atoms = len(atoms)
            
            print(f"  + Packed {n_waters} waters in {elapsed:.2f}s")
            print(f"  + Total atoms: {total_atoms} ({n_substrate} substrate + {n_waters*3} water)")
            
            # Store for output
            atoms.info['system'] = system_name
            atoms.info['size'] = size
            atoms.info['n_waters'] = n_waters
            atoms.info['time_s'] = f"{elapsed:.2f}"
            all_structures.append(atoms)
            
            results.append({
                'system': system_name,
                'size': size,
                'n_substrate': n_substrate,
                'n_waters': n_waters,
                'time': elapsed,
                'total_atoms': total_atoms
            })
    
    # Write all structures to single extxyz file
    output_file = "mb_benchmark_results.xyz"
    print(f"\n{'='*60}")
    print(f"Writing all structures to {output_file}...")
    write(output_file, all_structures, format='extxyz')
    print(f"+ Saved {len(all_structures)} structures")
    
    # Summary table
    print(f"\n{'='*60}")
    print("Summary:")
    print(f"{'System':<12} {'Size':<8} {'Substrate':<12} {'Waters':<10} {'Time(s)':<10}")
    print("-"*60)
    for r in results:
        print(f"{r['system']:<12} {r['size']:<8} {r['n_substrate']:<12} {r['n_waters']:<10} {r['time']:<10.2f}")
    
    print(f"\n{'='*60}")
    print(f"++ Benchmark complete! Inspect structures with:")
    print(f"   ase gui {output_file}")
    print("="*60)


if __name__ == "__main__":
    run_benchmark()
