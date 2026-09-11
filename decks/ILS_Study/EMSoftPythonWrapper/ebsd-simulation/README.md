# EMsoft python handler
## Description
- run the test_parser.py skript and provide the path to the cif file and a parameter file
    - these files should be located within your em-soft work folder
    - the skript will automatically read the em-soft work folder base on the emsoftconfig.json path in the config.yml file

- the skript will automatically
    - generate a emsoft .xtal file based on the cif input (please report bugs, as cif is only mediocre standardise)
    - generate the required .bash and .nml files (for the bash files we still need to find a solution, to hinder that when pushing / pulling you overwrite the basic        configuration)
    - run the EMMCOpenCL and EMEBSDmaster as a bash script automatically
    - the git hash of the current commit will be written to the .txt file to retrieve changes in the software at a later change
    - TODO: load up the results and respecitve metadata to openbis

## Visuals
Depending on what you are making, it can be a good idea to include screenshots or even a video (you'll frequently see GIFs rather than actual videos). Tools like ttygif can help, but check out Asciinema for a more sophisticated method.

## Installation
```bash
git clone https://git.rwth-aachen.de/imm_git/python_template
cd python_template
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```
## Configuration

open config.yml
change the path to your local EMSfot .json file this file wile contain the path necessary for EMsoft as specified in https://github.com/EMsoft-org/EMsoft/wiki/Package-Configuration


## Usage
This code is supposed to be used as a command line tool:
run
```bash
python3 test_parser.py -file= path/to/your/parameter_file.yaml -cif =path/to/you/cif_file.cif
```
If you have setup your EMsoft installation correctly (espeically the EMSoft.json file in your .config folder). The script will automatically start the EMmkxtal, EMMCOpenCl and EMEBSDmaster programms. (via a bash process)
This worked on the RWTH HPC Cluster on Claix-18 and Claix-23

Please make sure, to addjust your jocbscript_mc.sh and jobscript_master.sh files in the support folder according to your needs. i. e. check how you install EMSoft and change the export PATH=~/EMsoftBuild/Release/Bin:$PATH statement.

Please also change the parameters in default_config_emsoft.yml (especiall under bashscripts).
e.g. chhange the e-mail to your personal email. the other parameters can be found from SLURM documentation and from the RWTH hpc documentation

## Roadmap
Inclue Open Bis support

## Authors and acknowledgment


## License
For open source projects, say how it is licensed. Adjust the license file.

## Project status
This project is still under development.
    -> I.e. for some space groups the origin is hard coded. Right now we try to issue a warning if done so.
