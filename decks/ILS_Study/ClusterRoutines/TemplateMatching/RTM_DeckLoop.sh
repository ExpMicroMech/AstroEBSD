#!/usr/bin/zsh

### Job name
#SBATCH --job-name=RTMOPtimization

### File / path which STDOUT will be written to, %J is the job ID
#SBATCH --output=multiRTMCopies/logs/MATLAB_SERIAL.%J
#SBATCH -e multiRTMCopies/logs/MATLAB_ERROR_%j.err

### Request the time you need for execution. The full format is D-HH:MM:SS
### You must at least specify minutes or days and hours and may add or
### leave out any other parameters
#SBATCH --nodes=1
#SBATCH --time=23:59:00
##SBATCH --time=00:59:00
#SBATCH --mem=128000M
##SBATCH --mem=20000M
##SBATCH --cpus-per-task=1
#SBATCH --account=rwth1308
##SBATCH --partition=devel
#SBATCH --ntasks=25
##SBATCH --ntasks=4
### Change to the work directory, if not in submit directory
#cd "/work/mz071159/astro_skripts/decks/PCA"

### Load the required modules
# module load MATLAB/2023b



## WARNING: SOFAR WE ONLY SPECIFY THE ROTATION MATRICES DIRECTLY IN THE MATLAB SCRIPT.



##
# Define the original MATLAB script
original_script="RTM_DECK_individual_Template.sh"
scriptFolder="multiRTMCopies"


# Define replacement strings
## slashes need to be escaped !!
## Folder setups

# ------------------------------------------------------------------------------------- ##
## Nb Ni
# phaseSubFolder=("NbNiLPR_58Nb42Ni" "NbNiLPR_51Nb_49Ni")
phaseSubFolder=("NbNi3aILGrid_LPR58Nb_42Ni_2026_01" "NbNi3aILGrid_LPR51Nb_49Ni_2026_01")
# phaseSubFolder=("NbNi3aILGrid_LPR_58Nb_42Ni_2026_02_OCC" "NbNi3aILGrid_LPR_51Nb_49Ni_2026_02_OCC")
# phaseSubFolder=("NbNiLPR_58Nb_42Ni_HighBethe" "NbNiLPR_51Nb_49Ni_HighBethe")
###### phaseSubFolder=("NbNi_SiteLatticeSearch_LPR58Nb_42Ni_2026_01" "NbNi_SiteLatticeSearch_LPR51Nb_49Ni_2026_01")
SampleFolder=("58Nb-42Ni_v2_septemberbatch" "51Nb-49Ni_septemberbatch")
pcRefineH5File=("NbNi_58Nb-Ni_V2_september_Area 6_OIM Map 1_PCRefinement___JOBID__gradients" "NbNi_Nb-49Ni_V2_september_Area 8_EBSD+EDS Map 1_PCRefinement___JOBID__gradients")
pcRefineOSCFile=("NbNi_58Nb-Ni_V2_september_Area 6_OIM Map 1" "NbNi_Nb-49Ni_V2_september_Area 8_EBSD+EDS Map 1")
EBSDFolder=("Area7" "Area10") 
pcrefineFolder=("Area6" "NbNi_Nb-49Ni_V2_september_Area 8_EBSD+EDS_Map_1")

#--------------------------------------------------------------------------------------------##
## Nb Co 
# phaseSubFolder=("NbCo3aILGrid_LPR_High_2026_01" "NbCo3aILGrid_LPR_Low_2026_01"  "NbCo3aILGrid_LPR_Low_2026_01" )
############## phaseSubFolder=("NbCo_SiteLatticeSearch_HighLPR_2026_01" "NbCo_SiteLatticeSearch_LowLPR_2026_01"  "NbCo_SiteLatticeSearch_LowLPR_2026_01" )
# phaseSubFolder=("NbCo3aILGrid_LPR_High_2026_02_OCC" "NbCo3aILGrid_LPR_Low_2026_02_OCC" "NbCo3aILGrid_LPR_Low_2026_02_OCC")
# phaseSubFolder=("NbCo3aILGrid_LPR_High_2026_02_OCC_HighBethe" "NbCo3aILGrid_LPR_Low_2026_02_OCC_HighBethe" "NbCo3aILGrid_LPR_Low_2026_02_OCC_HighBethe")

# SampleFolder=("NbCoUBC_56Nb44Co_vertical" "NbCoUBC_48Nb52Co_repolishOPS_newCut_vert" "49Nb-51Co_WeiLuo_HT2")
# pcRefineH5File=("NbCoUBC_56Nb44Co_vertical_Area 2_EBSD+EDS Map 1_PCRefinement___JOBID__gradients" "NbCoUBC_48Nb52Co_repolishOPS_newCut_vert_Area 1_EBSD+EDS Map 2_PCRefinement___JOBID__gradients" "NbCoUBC_49Nb51Co_WeiLuo_HT2_Area 1_EBSD+EDS Map 2_PCRefinement___JOBID__gradients")
# pcRefineOSCFile=("NbCoUBC_56Nb44Co_vertical_Area 2_EBSD+EDS Map 1" "NbCoUBC_48Nb52Co_repolishOPS_newCut_vert_Area 1_EBSD+EDS Map 2" "NbCoUBC_49Nb51Co_WeiLuo_HT2_Area 1_EBSD+EDS Map 2")
# pcrefineFolder=("Area2" "Area1" "Area1")
# EBSDFolder=("Area3" "Area2" "Area2") 
# 

### ----------------------------------------------- #####
## Y Ni
# phaseSubFolder=("YNi2_30kVVoltageVar" "YNi2_30kVVoltageVar") #  "NbCo3aILGrid_LPR_Low_2026_01" )
# ############## phaseSubFolder=("NbCo_SiteLatticeSearch_HighLPR_2026_01" "NbCo_SiteLatticeSearch_LowLPR_2026_01"  "NbCo_SiteLatticeSearch_LowLPR_2026_01" )
# # phaseSubFolder=("NbCo3aILGrid_LPR_High_2026_02_OCC" "NbCo3aILGrid_LPR_Low_2026_02_OCC" "NbCo3aILGrid_LPR_Low_2026_02_OCC")
# # phaseSubFolder=("NbCo3aILGrid_LPR_High_2026_02_OCC_HighBethe" "NbCo3aILGrid_LPR_Low_2026_02_OCC_HighBethe" "NbCo3aILGrid_LPR_Low_2026_02_OCC_HighBethe")

# SampleFolder=("YNi2_NS_D" "YNi2_S_D") # "49Nb-51Co_WeiLuo_HT2")
# pcRefineH5File=("MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3_PCRefinement___JOBID__gradients" "MA_Dinah_YNi2_S_D_repolish_Area 5_EBSD+EDS Map 3_PCRefinement___JOBID__gradients") #"NbCoUBC_49Nb51Co_WeiLuo_HT2_Area 1_EBSD+EDS Map 2_PCRefinement___JOBID__gradients")
# pcRefineOSCFile=("MA_Dinah_YNi2_NS_D_repolish_Area 4_EBSD+EDS Map 3" "MA_Dinah_YNi2_S_D_repolish_Area 5_EBSD+EDS Map 3") # "NbCoUBC_49Nb51Co_WeiLuo_HT2_Area 1_EBSD+EDS Map 2")
# pcrefineFolder=("Area4" "Area5") # "Area1")
# EBSDFolder=("Area4" "Area5" ) #  "Area2") 


#  #copy the pha files you want to use into a folder in your Astro EBSD phases/phasesfiles/(newFolder) location
# # phaseSubFolder="LaNi"
# ## make sure to escape slashes use syntax: \/
rootFolder="/hpcwork/mz071159/EBSD/" 
Vtfolder="vt0,1" # for now we want to use the same for all.

# Define the new MATLAB script name
rotationMatrixPseudo="rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(60*degree)" # for trigonal pseudosymmetry
# rotationMatrixPseudo="rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(0*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(0*degree),rtmConfigRefine.RTM.Rz(90*degree)*rtmConfigRefine.RTM.Rx(90*degree)*rtmConfigRefine.RTM.Rz(90*degree)" ##for cubic - tetragonal pseudosymmetry
# rotationMatrixPseudo=""
# new_script="RTM_PCA_CommercialHough_map${SLURM_JOB_ID}.m"
# script_name="RTM_PCA_CommercialHough_map${SLURM_JOB_ID}"



# Copy the original script to a new one

# MultiPCABashCopies




for (( i = 1; i <= 2; i++ )); do  # zsh: start with 1 !!!!!!

    # Use sed to replace placeholders with specified strings
    # Folder setups

    DIRECTORY="${rootFolder}/${SampleFolder[i]}/${EBSDFolder[i]}"
    echo "${SampleFolder[i]}"
    PATTERN="Map*"
    pushd $DIRECTORY
    echo "Displaying Directory \n"
    echo "$DIRECTORY"
    echo "$new_script"
    # Save the find command output into a variable
    RESULT=$(find .  -type d -maxdepth 1 -name "*$PATTERN*" -printf '%P\n')
    popd
    # # Check if any results were found
    if [ -z "$RESULT" ]; then
        echo "No files or directories found matching pattern '$PATTERN' in '$DIRECTORY'."
    else
        echo "Found files and directories:"
        echo "$RESULT"
    fi
    for folder in $RESULT; do
        DIRECTORY="${rootFolder}/${SampleFolder[i]}/${EBSDFolder[i]}/${folder}"
        pushd  $DIRECTORY
        EBSDMap=$(find .  -maxdepth 1 -name "*.ang" -execdir basename {} .ang ';')
        popd
        new_script="RTM_PCA_CommercialHough_map${phaseSubFolder[i]}_${SampleFolder[i]}_${folder}.sh"
        cp "$original_script" "$new_script"
        newFolder="multiRTMCopies"
        # echo "$file"
        # pcrefineFolderIn="${SampleFolder[i]}\/${pcrefineFolder[i]}"
        # PatternFolder="${SampleFolder[i]}\/${pcrefineFolder[i]}\\/${folder}"
        # echo "$PatternFolder"

        sed -i "s/__PHASESSUBFOLDER__/$phaseSubFolder[i]/g" "$new_script"
        sed -i "s/__EBSDFOLDERIN__/$folder/g" "$new_script"
        sed -i "s/__EBSDMAPIN__/$EBSDMap/g" "$new_script"
        # sed -i "s/__AREAFo/$folder/g" "$new_script"
        sed -i "s/__MAPFOLDER__/$folder/g" "$new_script"
        sed -i "s/__SAMPLEFOLDER__/$SampleFolder[i]/g" "$new_script"
        sed -i "s/__AREAFOLDER__/$EBSDFolder[i]/g" "$new_script"
        sed -i "s/__VTFOLDER__/$Vtfolder/g" "$new_script"
        # sed -i "s/__PatternFolder__/$PatternFolder/g" "$new_script"

        sed -i "s/__PCREFINEOSCIN__/$pcRefineOSCFile[i]/g" "$new_script"
        sed -i "s/__PCREFINEFILEIN__/$pcRefineH5File[i]/g" "$new_script"
        sed -i "s/__PCREFINEFOLDERIN__/$pcrefineFolder[i]/g" "$new_script"
        sbatch "$new_script"
        done
    done