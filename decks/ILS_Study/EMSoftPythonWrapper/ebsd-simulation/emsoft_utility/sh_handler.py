import yaml
from pathlib import *

# modification and setup of nml files


def modify_sh_master(filename, simfolder, params):
    """change sh file based on input parameters
    might be redundant with definitions from
    below change in the future"""

    with open("support/jobscrit_master.sh", "r") as f:
        shell_script = f.read()
    with open("support/default_config_emsoft.yml") as f:
        conf = yaml.safe_load(f)
    for key in conf["bashscripts"].keys():
        if key in params.keys():
            shell_script = shell_script.replace(f"__{key}__", str(params[key]))
        else:
            shell_script = shell_script.replace(
                f"__{key}__", str(conf["bashscripts"][key])
            )
    shell_script = shell_script.replace("__jobname__", f"{filename.stem}_master")
    shell_script = shell_script.replace(
        "__filepath__", f"{simfolder.joinpath(filename.stem)}_master_"
    )
    shell_script = shell_script.replace(
        "__nml_file__", f"{simfolder.joinpath(filename.stem)}_master.nml"
    )
    with open(simfolder.joinpath(f"{filename.stem}_master.sh"), "w") as f:
        f.write(shell_script)


def modify_sh_mc(filename, simfolder, params):
    """change sh file based on input parameters
    might be redundant with definitions from
    below change in the future"""

    with open("support/jobscript_mc.sh", "r") as f:
        shell_script = f.read()
    with open("support/default_config_emsoft.yml") as f:
        conf = yaml.safe_load(f)
    for key in conf["bashscripts"].keys():
        # print(key)
        # print(key in params.keys())
        if key in params.keys():
            shell_script = shell_script.replace(f"__{key}__", str(params[key]))
        else:
            shell_script = shell_script.replace(
                f"__{key}__", str(conf["bashscripts"][key])
            )
    shell_script = shell_script.replace("__jobname__", f"{filename.stem}_mc")
    shell_script = shell_script.replace(
        "__filepath__", f"{simfolder.joinpath(filename.stem)}_mc_"
    )
    shell_script = shell_script.replace(
        "__nml_file__", f"{simfolder.joinpath(filename.stem)}_mc.nml"
    )
    with open(simfolder.joinpath(f"{filename.stem}_mc.sh"), "w") as f:
        f.write(shell_script)


def modify_bash_run(filename, simfolder, params):
    with open("support/multibash.sh", "r") as f:
        shell_script = f.read()
    shell_script = shell_script.replace(
        "__mc_sh__", f"{simfolder.joinpath(filename.stem)}_mc.sh"
    )
    shell_script = shell_script.replace(
        "__master_sh__", f"{simfolder.joinpath(filename.stem)}_master.sh"
    )
    with open(simfolder.joinpath(f"{filename.stem}_multi.sh"), "w") as f:
        f.write(shell_script)
