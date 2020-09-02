import ase.spacegroup
from ase.calculators.lammps import Prism  # , convert
import ase.io
import sys
import os
import numpy as np
from .crystals import crystals
from .geometry import BoxGeometry, PlaneBoundTriclinicGeometry
import requests
import requests_cache
import tempfile
from shutil import copyfile
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

    # cellpar = [a, b, c, alpha, beta, gamma]
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
    xlo, ylo, zlo = 0, 0, 0
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


def fetch_prepared_system(name, type_mapping=None):
    """Retrieves molecular system from an online repository. Caches data locally with requests-cache, which creates a sqlite database locally.

    Arguments:
    Name -- name of system to be retrieved
    type_mapping -- List of pairs of numbers mapping the atom type in the file from the online repository to an atom number. For example silica from an online repository of lammps data files will have Si and O as type 1 and 2, whereas the correct atomic numbers are 14 and 8

    Returns
        atoms -- ase.Atoms object with the system
    """
    requests_cache.install_cache('python_molecular_builder_cache')
    f = tempfile.TemporaryFile(mode="w+t")
    url = f"https://zenodo.org/record/3994120/files/{secure_filename(name)}.data"
    print("URI: ", url)
    r = requests.get(url, stream=True)
    print(r.encoding)
    r.encoding = "utf-8"
    total_length = int(r.headers.get('content-length'))
    print(f"Downloading data file {name}")
    chunk_size = 4096
    for chunk in progress.bar(r.iter_content(chunk_size=chunk_size, decode_unicode=True), expected_size=(total_length/chunk_size) + 1):
        if chunk:
            f.write(chunk)
            f.flush()
    f.seek(0)

    atoms = ase.io.read(f, format="lammps-data", style="atomic")
    if type_mapping is not None:
        type_copy = atoms.get_atomic_numbers()
        for pair in type_mapping:
            type_copy[atoms.numbers == pair[0]] = pair[1]
        atoms.set_atomic_numbers(type_copy)
    return atoms


def pack_water(atoms=None, nummol=None, volume=None, density=0.997,
               geometry=None, side='in', pbc=0.0, tolerance=2.0):
    """Pack water molecules into voids at a given volume defined by a geometry.
    The packing is performed by packmol.

    :param atoms: ase Atoms object that specifies particles that water is to be packed around. The packed water molecules will be added to this atoms object.
    :type atoms: Atoms object
    :param nummol: Number of water molecules
    :type nummol: int
    :param volume: Void volume in :math:`Å^3` to be filled with water. Can only be used if `nummol=None`, since `pack_water` will compute the number of atoms based on the volume and density of water.
    :type volume: float
    :param density: Water density. Used to compute the number of water molecules to be packed if `volume` is provided.
    :type density: float
    :param geometry: Geometry object specifying where to pack water
    :type geometry: Geometry object
    :param side: Pack water inside/outside of geometry
    :type side: str
    :param pbc: Inner margin to add to the simulation box to avoid overlapping atoms over periodic boundary conditions. This is necessary because packmol doesn't support periodic boundary conditions.
    :type pbc: float or array_like
    :param tolerance: Minimum separation distance between molecules.
    :type tolerance: float

    :returns: Coordinates of the packed water
    """
    if (volume is None and nummol is None):
        raise ValueError("Either volume or the number of molecules needed")
    elif (volume is not None) and (nummol is not None):
        raise ValueError("Either volume or the number of molecules needed")

    if volume is not None:
        V_per_water = 29.9796/density
        nummol = int(volume/V_per_water)

    format_s, format_v = "pdb", "proteindatabank"
    side += "side"

    if atoms is None and geometry is None:
        raise ValueError("Either atoms or geometry has to be given")
    elif geometry is None:
        # The default water geometry is a box which capsules the solid
        if type(pbc) is list or type(pbc) is tuple:
            pbc = np.array(pbc)

        if atoms.cell.orthorhombic:
            ll_corner = np.array([0,0,0])
            ur_corner = atoms.cell.lengths()
            box_center = (ur_corner + ll_corner) / 2
            box_length = ur_corner - ll_corner
            geometry = BoxGeometry(box_center, box_length - pbc)
        else:
            geometry = PlaneBoundTriclinicGeometry(atoms.cell, pbc=pbc)
    else:
        ll_corner = geometry.ll_corner
        ur_corner = geometry.ur_corner
        box_center = (ur_corner + ll_corner) / 2
        box_length = ur_corner - ll_corner

    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)
        sys.path.append(tmp_dir)

        if atoms is not None:
            # Write solid structure to pdb-file
            atoms.write(f"atoms.{format_s}", format=format_v)

        # Copy water.pdb to templorary directory
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
                f.write(f"  fixed 0 0 0 0 0 0\n")
                f.write("end structure\n\n")
            f.write(geometry.packmol_structure(nummol, side))

        # Run packmol input script
        try:
            os.system("packmol < input.inp")
        except:
            raise OSError("packmol is not found. For installation instructions, \
                           see http://m3g.iqm.unicamp.br/packmol/download.shtml.")

        # Read packmol outfile
        water = ase.io.read(f"out.{format_s}", format=format_v)

    os.chdir(cwd)

    if atoms is not None:
        # Remove solid
        del water[:len(atoms)]

    # Scale water box correctly
    #print(water.cell)
    #print(np.diag(box_length))
    #water.set_cell(np.diag(box_length))
    if atoms is not None:
        atoms += water

    return water


def write(atoms, filename, bond_specs=None, atom_style="molecular", size=(640, 480), camera_dir=(2, 1, -1), viewport_type="perspective", atom_radii=None):
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

        types = pipeline.source.data.particles.particle_types
        for symbol, i in symbols_dict.items():
            types.type_by_id(i).name = symbol
            types.type_by_id(i).load_defaults()


        if not atom_radii is None:
            types = pipeline.source.data.particles.particle_types
            for pair in atom_radii:
                types.type_by_id(symbols_dict[pair[0]]).radius = pair[1]

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
            export_file(pipeline, filename, "lammps/data", atom_style=atom_style)

        elif suffix == ".png":

            from ovito.vis import Viewport, TachyonRenderer, OpenGLRenderer
            pipeline.add_to_scene()

            if viewport_type =="perspective":
                vp = Viewport(type = Viewport.Type.Perspective, camera_dir = camera_dir)
            elif viewport_type =="orthogonal":
                vp = Viewport(type = Viewport.Type.Ortho, camera_dir = camera_dir)
            else:
                raise ValueError("viewport type has to be perspective or orthogonal")

            vp.zoom_all(size=size)
            vp.render_image(filename=filename, size=size, renderer=TachyonRenderer())
            pipeline.remove_from_scene()
