# import numpy as np
# from gemmi import cif
# import pandas as pd
from emsoft_utility import cif_handler, debye_waller, nml_handler, sh_handler
from pathlib import *
import subprocess
import shutil

# import sys
import argparse
import yaml

# import os
import re


def gen_emsoft_txt(filename, folder, dbw=None):
    # print(f"dbw: {dbw}")
    lattice, atoms = cif_handler.read_cif(filename)
    atoms_debye = debye_waller.assign_debye(atoms, dbw)
    lattice["cs_emsoft"] = cif_handler.assign_crystal_system(lattice["space_group"])
    Path("dump").mkdir(parents=True, exist_ok=True)
    cif_handler.write_head(f"{folder.joinpath(filename.stem)}.txt", lattice)
    cif_handler.write_foot(f"{folder.joinpath(filename.stem)}.txt", atoms_debye)
    return [lattice, atoms]


def main():
    #    filename='fe6ta7.cif'
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-cif",
        "--ciffile",
        help="path of the cif file respective to your em-soft work directory from which you want to start the simulation",
        type=str,
        default="setup/fe6ta7/fe6ta7.cif",
    )
    parser.add_argument(
        "-pfile",
        "--paramfile",
        help="path to parameter file where you change with respect to default parameters (relative to em-soft work directory)",
        type=str,
    )
    parser.add_argument(
        "-dbw",
        "--overwrite_waller",
        help='overwrite the debeye waller parameters -> might lead to "crisper" patterns, emsoft recommends [0.004 - 0.006]',
        type=float,
        # default=0.005,
    )
    args = parser.parse_args()

    filename = args.ciffile
    filename = Path(nml_handler.get_emdatapath()).joinpath(filename)
    parfile = args.paramfile
    if not parfile is None:
        parfile = Path(nml_handler.get_emdatapath()).joinpath(parfile)
        print(parfile)
        print(nml_handler.get_emdatapath())
        with open(parfile, "r") as f:
            params = yaml.safe_load(f)
    else:
        params = {}
    print(params)
    simfolder = nml_handler.setup_nml_files_default(filename, params)
    # gernerate the xtal files
    # print(f"dbwalleroverwrite:{args.overwrite_waller}")
    if args.overwrite_waller is None:
        [lattice, atoms] = gen_emsoft_txt(Path(filename), simfolder)
    else:
        [lattice, atoms] = gen_emsoft_txt(
            Path(filename), simfolder, args.overwrite_waller
        )
    input = open(f"{simfolder.joinpath(filename.stem)}.txt")
    subprocess.call("EMmkxtal", stdin=input)
    # a
    # setup simulation

    # print(simfolder)
    sh_handler.modify_sh_mc(filename, simfolder, params["sh"])
    sh_handler.modify_sh_master(filename, simfolder, params["sh"])
    print(params)

    # generate .ph files
    phaseFileName = f"{filename.stem}_{simfolder.parents[0].stem}_{simfolder.stem}.pha"
    cif_handler.writePha(simfolder.joinpath(phaseFileName), lattice, filename, params)

    ## copy cif files and create .pha file
    print(filename)
    print(simfolder.joinpath(f"{filename.stem}.cif"))
    shutil.copy(filename, str(simfolder.joinpath(f"{filename.stem}.cif")))
    # ##run simulation
    cmd = ["sbatch", f"{simfolder.joinpath(filename.stem)}_mc.sh"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    cmd_out, cmd_err = p.communicate()
    print(cmd_out)
    m = re.search("(?<=job )\d*", str(cmd_out))
    if not m is None:
        jobid = m.group(0)
        cmd = [
            "sbatch",
            f"--dependency=afterok:{jobid}",
            f"{simfolder.joinpath(filename.stem)}_master.sh",
        ]
        subprocess.call(cmd)


if __name__ == "__main__":
    print("running main")
    main()
