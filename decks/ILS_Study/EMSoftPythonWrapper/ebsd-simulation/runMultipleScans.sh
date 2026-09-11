####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.250_c26.7_a4.9.cif
####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.275_c26.7_a4.9.cif
####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.300_c26.7_a4.9.cif
####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.325_c26.7_a4.9.cif
####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.350_c26.7_a4.9.cif
####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.375_c26.7_a4.9.cif
####python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.400_c26.7_a4.9.cif
########## python test_parser.py -pfile=setup/NbNi_mu/params.yml -cif=setup/NbNi_mu/Nb553Ni447_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_18.0kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_18.5kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_18.8kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_19.0kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_19.8kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_20.0kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_20.2kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_20.5kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_21.0kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
###python test_parser.py -pfile=setup/NbNi_mu/params_20.2kV.yml -cif=setup/NbNi_mu/Nb518Ni482_joubert.cif
##python test_parser.py -pfile=setup/RENi2-Laves/params.yml -cif=setup/RENi2-Laves/Y95Ni2.cif
# python test_parser.py -pfile=setup/RENi2-Laves/Y/30kV.yml -cif=setup/RENi2-Laves/Y/YNi2.cif
##python test_parser.py -pfile=setup/RENi2-Laves/params.yml -cif=setup/RENi2-Laves/ErNi2.cif
##python test_parser.py -pfile=setup/RENi2-Laves/params.yml -cif=setup/RENi2-Laves/Er60Ni132.cif
##python test_parser.py -pfile=setup/RENi2-Laves/params.yml -cif=setup/RENi2-Laves/YNi2.cif
###python test_parser.py -pfile=setup/Meteorite_SangHyeok/setup.yml -cif=setup/Meteorite_SangHyeok/tetrataenite_vesta.cif

#### python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.425_c26.7_a4.9.cif
#### python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.450_c26.7_a4.9.cif
#### python test_parser.py -pfile=setup/NbNi_mu/params_19.5kV.yml -cif=setup/NbNi_mu/NbNi_dt0.475_c26.7_a4.9.cif


#!/bin/bash

# Define the path to your parameter file
# PARAM_FILE="/hpcwork/mz071159/EMsoft_work/setup/NbNi_mu/3aMixedEnergySeries"
# PARAM_FILE="/hpcwork/mz071159/EMsoft_work/setup/Meteorite_SangHyeok/setup.yml"
# PARAM_FILE="/hpcwork/mz071159/EMsoft_work/setup/NbNi_mu/DMinVar/params_19.5kV_0.05.yml"
# PARAM_FILE="/hpcwork/mz071159/EMsoft_work/setup/Ni/params_19.5kV.yml" # /params_19.5kV_HighBethe.yml" #/NbCo_mu/params_19.5kV.yml"
PARAM_FILE="/hpcwork/mz071159/EMsoft_work/setup/RENi2-Laves/Y/YNi_param_var"
# PARAM_FILE="/hpcwork/mz071159/EMsoft_work/setup/stacking_faults/params.yml" # /params_19.5kV_HighBethe.yml" #/NbCo_mu/params_19.5kV.yml"
 # Define the directory containing .cif filer
# CIF_DIR="/hpcwork/mz071159/EMsoft_work/setup/NbNi_mu/LPRMeasuredNew_2026_01_16_IL"

# CIF_DIR="/hpcwork/mz071159/EMsoft_work/setup/NbNi_mu/DMinVar"
CIF_DIR="/hpcwork/mz071159/EMsoft_work/setup/Ni" #/LPRMeasuredNew_2026_03_05_IL"
# CIF_DIR="/hpcwork/mz071159/EMsoft_work/setup/NbNi_mu/LPRMeasuredNew_2026_02_09_OCC"
CIF_DIR="/hpcwork/mz071159/EMsoft_work/setup/stacking_faults" #/LPRMeasuredNew_2026_03_05_IL"
# CIF_DIR="/hpcwork/mz071159/EMsoft_work/setup/Si"
CIF_FILE="setup/RENi2-Laves/Y/Y95Ni2.cif"
# Loop through all .cif files in the specified directory
# for cif_file in "$CIF_DIR"/*.cif; do
#     if [ -f "$cif_file" ]; then  # Check if it is a file
#         echo "Running command for: $cif_file"
#         python test_parser.py -pfile="$PARAM_FILE" -cif="$cif_file"
#     fi
# done

# Loop through all .yml files in the specified directory
for par_file in "$PARAM_FILE"/*.yml; do
    if [ -f "$par_file" ]; then  # Check if it is a file
        echo "Running command for: $par_file"
        python test_parser.py -pfile="$par_file" -cif="$CIF_FILE" #-dbw=0.0002
    fi
done