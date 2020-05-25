import ase.spacegroup
from ase.calculators.lammps import Prism, convert
import ase.io
import sys
import os
import numpy as np
from .crystals import crystals 
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
    


def pack_water(atoms, number, geometry, side='inside', pbc=None, tolerance=2.0):
    """Pack water molecules into voids at a given volume defined by a geometry.
    
    :param atoms: ase Atom object that specifies where the solid is
    :type atoms: Atom object
    :param number: Number of water molecules
    :type number: int
    :param geometry: Geometry object specifying where to pack water
    :type geometry: Geometry object
    
    :returns: Coordinates of the packed water
    """
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        sys.path.append(tmp_dir)
        
        # Write solid structure to pdb-file
        atoms.write("atoms.pdb", format="proteindatabank")
        
        # ...
        from shutil import copyfile
        this_dir, this_filename = os.path.split(__file__)
        water_data = this_dir + "/data_files/water.pdb"
        copyfile(water_data, "water.pdb")
        
        # Generate packmol input script
        with open("input.inp", "w") as f:
            f.write(f"tolerance {tolerance}\n")
            f.write("filetype pdb\n")
            f.write("output out.pdb\n")
            f.write("structure atoms.pdb\n")
            f.write("  number 1\n")
            f.write("  center\n")
            f.write("  fixed 0 0 0 0 0 0\n")
            f.write("end structure\n\n")
            f.write(geometry.packmol_structure(number, side))
        
        # Run packmol input script
        try:
            os.system("packmol < input.inp")
        except:
            raise OSError("packmol is not found. For installation instructions, see http://m3g.iqm.unicamp.br/packmol/download.shtml.")
        
        # Read packmol outfile
        solid_and_water = ase.io.read("out.pdb", format="proteindatabank")
    return solid_and_water

