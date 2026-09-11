# Cif handling
from gemmi import cif
import numpy as np
import shutil

# from gemmi import cif
import pandas as pd
from pathlib import Path
import git


def read_cif(filename):
    """this reads in the lattice parameters and space group from the cif file, that we need to put in to the EM soft simulation
        Right now we assume, that the relevant information is in the first block of the cif file, which is at least true for vest files,
            if bugs occur with this setting please report and we will try to fix"""
    doc = cif.read_file(str(filename))
    block = doc.sole_block()
    # get the lattice parameters
    lattice = {}
    lattice["a"] = to_float(block.find_value("_cell_length_a"))
    lattice["b"] = to_float(block.find_value("_cell_length_b"))
    lattice["c"] = to_float(block.find_value("_cell_length_c"))
    lattice["alpha"] = to_float(block.find_value("_cell_angle_alpha"))
    lattice["beta"] = to_float(block.find_value("_cell_angle_beta"))
    lattice["gamma"] = to_float(block.find_value("_cell_angle_gamma"))
    lattice["space_group"] = int(block.find_value("_space_group_IT_number"))
    # atoms with seperate skript
    atoms = get_atom_positions(block)

    return lattice, atoms


                                                                    # def read_contcar(filename)


def get_atom_positions(block):
    """Get all the atom positions and add them to a pandas Dataframe"""
    atom_pos = {}
    atom_pos["element"] = list(block.find_loop("_atom_site_type_symbol"))
    atom_pos["x"] = list(block.find_loop("_atom_site_fract_x"))
    atom_pos["y"] = list(block.find_loop("_atom_site_fract_y"))
    atom_pos["z"] = list(block.find_loop("_atom_site_fract_z"))
    atom_pos["occ"] = list(block.find_loop("_atom_site_occupancy"))
    # atom_pos['occ']=np.ones(len(atom_pos['element']))''
    u_iso = list(block.find_loop("_atom_site_U_iso_or_equiv"))
    print(u_iso)
    bool_iso = not (u_iso)
    print(v == 0 for v in u_iso)
    if (
        not any([v == 0 for v in u_iso])
        and not bool_iso
        and all([v.isnumeric() for v in u_iso])
    ):
        atom_pos["dbw"] = u_iso
    print(atom_pos)
    atomic_positions = pd.DataFrame.from_dict(atom_pos)
    return atomic_positions

# utitlity
def to_float(input):
    """Convert the string inputs to int, if they contain a round bracket (in most cases denoting an error value in lattice parameter measurements)"""
    if not input.find("(") == -1:
        input = input[0 : input.find("(")]
        input = float(input)
        return input
    else:
        return float(input)


def assign_crystal_system(space_group):
    """Return the crystal system number based on the EMmkxtal command and comparison with internation space group table number
    1. Cubic
    2. Tetragonal
    3. Orthorhombic
    4. Hexagonal
    5. Trigonal
    6. Monoclinic
    7. Triclinic

    """
    space_dict = {
        7: (1, 2),
        6: (3, 15),
        3: (16, 74),
        2: (75, 142),
        5: (143, 167),
        4: (168, 194),
        1: (195, 230),
    }
    for system, space_group_val in space_dict.items():
        if space_group_val[0] <= space_group <= space_group_val[-1]:
            return system


def write_head(filename, lattice):
    """writes the necessary inputs for the txt file for the emsoft to generate the .xtal file"""

    with open(filename, "w") as f:
        sg = lattice["space_group"]
        cs = lattice["cs_emsoft"]
        a = lattice["a"] / 10  # em soft uses nanometers
        b = lattice["b"] / 10
        c = lattice["c"] / 10
        alpha = lattice["alpha"]
        beta = lattice["beta"]
        gamma = lattice["gamma"]

        f.write(f"{cs}\n")
        if cs == 1:  # cubic
            f.write(f"{a}\n")
        elif cs == 2:  # tetragonal
            f.write(f"{a}\n")
            f.write(f"{c}\n")
        elif cs in [4, 5]:  # hexagonal and trigonal
            if cs == 5:
                f.write("0\n")  # 0= non rhomboedal
            f.write(f"{a}\n")
            f.write(f"{c}\n")
        elif cs in [3, 6]:  # orthoromibc and monocliinic
            f.write(f"{a}\n")
            f.write(f"{b}\n")
            f.write(f"{c}\n")
            if cs == 6:  # angle for monolcinic
                f.write(f"{beta}\n")
        elif cs == 7:  # triclinic
            f.write(f"{a}\n")
            f.write(f"{b}\n")
            f.write(f"{c}\n")
            f.write(f"{alpha}\n")
            f.write(f"{beta}\n")
            f.write(f"{gamma}\n")
        f.write(f"{sg}\n")
        print(sg)
        if sg == 227:
            f.write("1\n")
            print(
                "warning, for this spacegroup we have hard coded the origin selection !!!!!!"
            )
        if sg == 59:
            f.write("1\n")
            print(
                "warning, for this spacegroup we have hard coded the origin selection !!!!!!"
            )


def write_foot(filename, atoms_debye):
    """writes the atom coordinates based on the input dataframe"""
    with open(filename, "a") as f:
        # print(atoms_debye)
        for i in atoms_debye.index:
            z = atoms_debye.loc[i, "Z"]
            f.write(f"{z}\n")
            x = atoms_debye.loc[i, "x"]
            y = atoms_debye.loc[i, "y"]
            z = atoms_debye.loc[i, "z"]
            occ = atoms_debye.loc[i, "occ"]
            dbw = atoms_debye.loc[i, "dbw"]
            f.write(f"{x},{y},{z},{occ},{dbw}\n")
            if i < len(atoms_debye.index) - 1:
                f.write("y\n")
            else:
                f.write("n\n")
        f.write(f"{Path(filename).stem}.xtal\n")

        repo = git.Repo(search_parent_directories=True)
        sha = repo.head.object.hexsha
        f.write(f"Source: undefined, git_hash: {sha}")


def writePha(filename, lattice, cifName,params):
    """writes the phase phile in AstroEBSD convention
    Our focus is to include the correct name of the masterfile, the cif file and the type of simulation. The rest will be generic, as it won't be use.
    We will include, whether the phase is hexagonal
    """
    cs = lattice["cs_emsoft"]
    with open(filename, "w") as f:
        # print(atoms_debye)
        f.write("$Name\n")
        f.write(f"{cifName.stem}\n")
        f.write("$Atoms\n")
        f.write("0.0000,0.0000,0.0000;0.5000,0.50000,0.50000\n")
        f.write("$Symmetry\n")
        f.write("2\n")
        f.write("$Lattice\n")
        f.write(f"{lattice['a']},{lattice['b']},{lattice['c']},{lattice['alpha']},{lattice['beta']},{lattice['alpha']},90\n")
        f.write("$FFactor\n")
        f.write("2,2,2,2\n")
        f.write("$MaxFam\n")
        f.write("8\n")
        f.write("$MaxIndex\n")
        f.write("20\n")
        f.write("$EFac\n")
        f.write("100\n")
        f.write("$PlaneReflector\n")
        f.write("0,0,1;1,0,0;1,0,3;0,1,3;1,1,1;1,1,3;\n")
        f.write("$cif\n")
        f.write(f"{cifName.stem}.cif\n")
        f.write("$sim_type\n")
        f.write("Emsoft\n")
        f.write("$dynamic\n")
        f.write(f"{filename.stem}.h5\n")
        f.write("$isHex\n")
        if cs in [4, 5]:  # hex and trigonal
            f.write("1\n")
        else:
            f.write("0\n")

        for key in params.keys():
            if not key == 'sh':
                for subkey in params[key].keys():
                    if subkey == 'ehistmin':
                        f.write('$Voltage\n')
                        f.write(f"{params[key][subkey]}\n")
                    f.write(f'${subkey}\n')
                    f.write(f"{params[key][subkey]}\n")

        f.write("$Comments\n")
        f.write(
            "WARNING: THIS FILE DOES NOT CONTAIN ACCURATE LATTICE INFORMATION. IT IS ONLY MEANT TO CONNET CIF FILE AND DYNAMIC TEMPLATE\n"
        )

    #  atoms = cif_handler.read_cif(filename)
