import ase.spacegroup
from ase.calculators.lammps import Prism, convert
import ase.io
import numpy as np
from .crystals import crystals 
import requests 
import requests_cache
import tempfile
from clint.textui import progress
from werkzeug.utils import secure_filename


def create_bulk_crystal(name, size, round="up"):
    """Create a bulk crystal from a spacegroup description.

    :param name: name of the crystal. A list can be found by @TODO
    :type name: str
    :param size: size of the bulk crystal. In the case of a triclinic cell, the dimensions are the ones along the diagonal of the cell matrix, and the crystal tilt decides the rest.  
    :type size: array_like with 3 elements

    :return: ase.Atoms object containing the crystal
    :rtype: ase.Atoms
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

def write(atoms, filename, bond_specs = None, size=(640, 480), atom_style="molecular"):
    """Write atoms to lammps data file 

    :param atoms: The atoms object to write to file 
    :type atoms: ase.Atoms 
    :param filename: filename to write to. Can either have suffix `.data` or `.png`, in which case a Lammps data file or a png picture will be produced, respectively. 
    :type filename: str 
    :param bonds_spec: List of (element1, element2, cutoff)
    :type bonds_spec: List of tuples
    """
    import os 
    suffix = os.path.splitext(filename)[1]

    if not suffix in [".data", ".png"]:
        raise ValueError(f"Invalid file format {suffix}")

    import tempfile
    import os
    from ase.formula import Formula 
    # Using tempfile and write + read rather than ovito's ase_to_ovito and back because the 
    # ordering of particle types for some (bug) reason becomes unpredictable 
    with tempfile.TemporaryDirectory() as tmp_dir:
        symbols = list(Formula(atoms.get_chemical_formula()).count().keys())
        symbols_dict = {}
        for i, symbol in enumerate(symbols): 
            symbols_dict[symbol] = i+1
        atoms.write(os.path.join(tmp_dir, "tmp.data"), format="lammps-data", specorder = symbols)
        
        from ovito.io import import_file, export_file
        from ovito.modifiers import CreateBondsModifier
        
        pipeline = import_file(os.path.join(tmp_dir, "tmp.data"))
        
        # Accept a single tuple not contained in a list if there is only one bond type. 
        if not bond_specs is None: 
            bondsmodifier = CreateBondsModifier(mode = CreateBondsModifier.Mode.Pairwise)
            if not isinstance(bond_specs, list) and isinstance(bond_specs, tuple):
                bond_specs = [bond_specs]
            for element in bond_specs:
                bondsmodifier.set_pairwise_cutoff(  symbols_dict[element[0]], 
                                                    symbols_dict[element[1]], 
                                                    element[2])
            pipeline.modifiers.append(bondsmodifier)
        pipeline.compute()
        if suffix == ".data":
            export_file(pipeline, filename, "lammps/data", atom_style=atoms_style)

        elif suffix == ".png":
            from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
            pipeline.add_to_scene()
            vp = Viewport(type = Viewport.Type.Perspective, camera_dir = (2, 1, -1))
            vp.zoom_all(size=size)
            vp.render_image(filename=filename, size=size, renderer=TachyonRenderer())
