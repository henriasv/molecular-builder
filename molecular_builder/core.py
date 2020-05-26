import ase.spacegroup
from ase.calculators.lammps import Prism, convert
import ase.io
import sys
import os
import numpy as np
from .crystals import crystals 
from .geometry import BoxGeometry
import requests 
import requests_cache
import tempfile
from clint.textui import progress
from werkzeug.utils import secure_filename


def create_bulk_crystal(name, size, round="up"):
    """Create a bulk crystal from a spacegroup description.

    Arguments: 
        name -- name of the crystal. A list can be found by @TODO
        size -- size of the bulk crystal. In the case of a triclinic cell, the dimensions are the ones along the diagonal of the cell matrix, and the crystal tilt decides the rest. 

    Returns:
        ase.Atoms object containing the crystal
    """
    crystal = crystals[name]
    a, b, c, alpha, beta, gamma = [crystal[i] for i in ["a", "b", "c", "alpha", "beta", "gamma"]]
    lx, ly, lz = size[0], size[1], size[2]
    
    cellpar = [a, b, c, alpha, beta, gamma]
    repeats = [lx/a, ly/b/np.sin(np.radians(gamma)), lz/c/np.sin(np.radians(alpha))/np.sin(np.radians(beta))]
    if round == "up":
        repeats = [int(np.ceil(i)) for i in repeats]
    elif round == "down":
        repeats = [int(np.floor(i)) for i in repeats]
    elif round == "round":
        repeats = [int(round(i)) for i in repeats]
    else:
        raise ValueError
    myCrystal = ase.spacegroup.crystal(
                    crystal["elements"],
                    crystal["positions"],
                    spacegroup = crystal["spacegroup"],
                    cellpar = [crystal[i] for i in ["a", "b", "c", "alpha", "beta", "gamma"]],
                    size=repeats)

    ###############################################################################
    # Creating a Lammps prism and then recreating the ase cell is necessary 
    # to avoid flipping of the simulation cell when outputing the lammps data file
    # By making the transformation here, what we see in the lammps output is the same as
    # the system we are actually carving into
    p = Prism(myCrystal.cell)
    xhi, yhi, zhi, xy, xz, yz = p.get_lammps_prism()
    xlo = 0; ylo = 0; zlo= 0
    cell = np.zeros((3, 3))
    cell[0, 0] = xhi - xlo
    cell[1, 1] = yhi - ylo
    cell[2, 2] = zhi - zlo
    if xy is not None:
        cell[1, 0] = xy
    if xz is not None:
        cell[2, 0] = xz
    if yz is not None:
        cell[2, 1] = yz

    myCrystal.set_cell(cell)
    myCrystal.wrap()
    ##################################################################################
    return myCrystal

def carve_geometry(atoms, geometry, side="in", return_carved=False):
    """Delete atoms according to geometry.

    Arguments: 
        atoms -- The ase Atoms object containing the molecular system 
        geometry -- A molecular_builder.geometry.Geometry object defining the region to be carved
        side -- Whether to carve out the inside or the outside of geometry

    Returns:  
        Number of deleted atoms 
        Optionally an atoms object containing the atoms that were carved away
    """

    if return_carved:
        atoms_copy = atoms.copy()
    
    geometry_indices = geometry(atoms)
    
    if side == "in":
        delete_indices = geometry_indices
    elif side == "out":
        delete_indices = np.logical_not(geometry_indices)
    else: 
        raise ValueError

    del atoms[delete_indices] 
    
    if not return_carved:
        return np.sum(delete_indices)
    else: 
        del atoms_copy[np.logical_not(delete_indices)]
        return np.sum(delete_indices), atoms_copy



def fetch_prepared_system(name):
    """Retrieves molecular system from an online repository. Caches data locally with requests-cache, which creates a sqlite database locally. 

    args: 
    Name -- name of system to be retrieved

    Returns 
        atoms -- ase.Atoms object with the system 
    """
    requests_cache.install_cache('python_molecular_builder_cache')
    f = tempfile.TemporaryFile(mode="w+t")
    url = f"https://zenodo.org/record/3774915/files/{secure_filename(name)}.data"
    print("URI: ", url)
    r = requests.get(url, stream=True)
    print(r.encoding)
    r.encoding="utf-8"
    total_length = int(r.headers.get('content-length'))
    print(f"Downloading data file {name}")
    chunk_size = 4096
    for chunk in progress.bar(r.iter_content(chunk_size=chunk_size, decode_unicode=True), expected_size=(total_length/chunk_size) + 1): 
        if chunk:
            f.write(chunk)
            f.flush()
    f.seek(0)

    atoms = ase.io.read(f, format="lammps-data", style="atomic")
    return atoms 
    


def pack_water(number, atoms=None, geometry=None, side='in', pbc=False, tolerance=2.0):
    """Pack water molecules into voids at a given volume defined by a geometry.
    
    :param number: Number of water molecules
    :type number: int
    :param atoms: ase Atoms object that specifies where the solid is
    :type atoms: Atoms object
    :param geometry: Geometry object specifying where to pack water
    :type geometry: Geometry object
    :param side: Pack water inside/outside of geometry
    :type side: str
    :param pbc: Ensures that the water molecules are separated by a certain distance through periodic boundaries. If float, the same distance is applied in all directions
    :type pbc: float or ndarray
    :param tolerance: minimum separation distance between molecules. 2.0 by default.
    :type tolerance: float
    
    :returns: Coordinates of the packed water
    """
    
    format_s, format_v = "pdb", "proteindatabank"    
    side += "side"
    
    # Geometrical properties of solid
    positions = atoms.get_positions()
    ll_corner = np.min(positions, axis=0)
    ur_corner = np.max(positions, axis=0)
    center = (ur_corner + ll_corner) / 2
    length = ur_corner - ll_corner
    
    if atoms is None and geometry is None:
        raise ValueError("Either atoms or geometry has to be given")
    elif geometry is None:
        # The default water geometry is a box which capsules the solid
        if pbc:
            center -= pbc / 2
            length -= pbc
        geometry = BoxGeometry(center, length)
    
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        sys.path.append(tmp_dir)
        
        if atoms is not None:
            # Write solid structure to pdb-file
            atoms.write(f"atoms.{format_s}", format=format_v)
        
        # Copy water.pdb to templorary directory
        from shutil import copyfile
        this_dir, this_filename = os.path.split(__file__)
        water_data = this_dir + f"/data_files/water.{format_s}"
        copyfile(water_data, f"water.{format_s}")
        
        # Generate packmol input script
        with open("input.inp", "w") as f:
            f.write(f"tolerance {tolerance}\n")
            f.write(f"filetype {format_s}\n")
            f.write(f"output out.{format_s}\n")
            if atoms is not None:
                f.write(f"structure atoms.{format_s}\n")
                f.write("  number 1\n")
                f.write("  center\n")
                f.write(f"  fixed {center[0]} {center[1]} {center[2]} 0 0 0\n")
                f.write("end structure\n\n")
            f.write(geometry.packmol_structure(number, side))
        
        # Run packmol input script
        try:
            os.system("packmol < input.inp")
        except:
            raise OSError("packmol is not found. For installation instructions, see http://m3g.iqm.unicamp.br/packmol/download.shtml.")
        
        # Read packmol outfile
        water = ase.io.read(f"out.{format_s}", format=format_v)
        
    os.chdir(cwd)
        
    # Remove solid
    del water[:len(atoms)]
    
    # Geometrical properties of water
    ll_corner = geometry.ll_corner
    ur_corner = geometry.ur_corner
    center = (ur_corner + ll_corner) / 2
    length = ur_corner - ll_corner
    
    # Scale water box correctly
    water.set_cell(np.identity(3) * length, scale_atoms=True)
    #water.positions *= length
    #water.positions += center
    return water

