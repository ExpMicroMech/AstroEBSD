AstroECP v1.1 - 2025

-- REFERENCE --
The AstroECP GUI is described in:
M. Haroon Qaiser, Lukas Berners, Robin J. Scales, Tianbi Zhang, Martin Heller, Jiří Dluhošd, Sandra Korte-Kerzel and T. Ben Britton, “AstroECP: Towards more practical Electron Channeling Contrast Imaging”, (2025) Preprint: arXiv (https://arxiv.org/abs/2507.00354)

-- ABOUT --
A digital twin of our SEM within a graphical user interface to help navigate the electron channeling patterns and quickly access specific ECCI conditions for subsequent analysis, based on the approaches and functions used in AstroEBSD (Britton, Tong et al., 2018) and MTEX v5.10.2 (Bachmann et al., 2010; Hielscher et al., 2019). In essence, this means that AstroECP serves not only as a pattern indexing and visualization tool but also as an experimental control assistant that facilitates ECCI by precisely dictating about the requisite stage tilts/rotations of SEM to navigate along the crystal.

-- HOW IT WORKS --
AstroECP reads the dynamically calculated Kikuchi reference patterns generated from pattern simulation software that include: MapSweeper EBSD dynamical patterns (higher quality Bloch Wave Kikuchi Diffraction (BWKD) patterns), EMsoft (Singh, Marc de Graefe et al., 2017) and Bruker DynamicS. The GUI then reprojects the channeling pattern based upon the crystal orientation and ECP conditions. Manipulation of the simulation can be performed through virtual movement of the stage and microscope to inform movement of the sample  within the SEM. This helps to effectively orient the crystal for the development of crystallographic contrast with a specific [uvw] direction along the optic axis of the SEM, and provide imaging conditions that are suitable for the collection of ECCI micrographs.

-- SETUP REQUIREMENTS --
To use AstroECP for pattern indexing, simulation and navigation, the following setup is required:
 - A MATLAB environment (R2020b or newer) with the MTEX toolbox (version 5.10.2 or later) and the AstroEBSD setup
 - A representative crystallographic information file (.cif) for the sample, containing essential unit cell parameters and symmetry data necessary for accurate pattern simulation and phase definition (an example of Si is included here).
 - A phase file (.pha) for the material of interest, which also specifies the relevant phase reflectors for kinematic band labelling and provides the path to simulated patterns. AstroEBSD includes several example phase files, and users may generate new ones following the provided template for different materials (an example of Si is included here).
 - A dynamically simulated reference pattern generated using BWKD (output as .h5 files), MapSweeper (as .sdf5 format), EMsoft (also .h5 format), or Bruker DynamicS (.bin files).
 - An experimental channeling pattern acquired from the SEM, typically in .tif format, accompanied by a text based header file (.hdr) that contains relevant microscope settings such as accelerating voltage, working distance, pattern size, stage positions etc. For instruments that do not create a text based header file, a user could create a similar file from one of the examples provided here. 

-- HOW TO USE IN MATLAB --
After setting up the files as mentioned in SETUP REQUIREMENTS, the user can open the AstroECP_v1.m and setup the paths and voltage:
Input_Data.astro_location='your AstroEBSD folder path';
Input_Data.mtex_location='your MTEX folder path'; %working with 5.10.2
Input_Data.image_frame=1; 
Input_Data.V_in = 20; 

Run AstroECP_v1.m as usual. This should give a full screen GUI.


-- FEATURES --
In addition to projection, initial crystal orientation, and subsequent stage tilt/rotations, the GUI also provides a simple unit cell visualization tool, as well as an overlay of the kinematic band edges to enable Miller family-based indexing of the bands within the Kikuchi pattern. There are a few other features within this software, including phase selection to compare different simulations, pattern matching based refinement of crystal orientation and/or projection parameters, enhanced ECP visualization (i.e. contrast stretching within the histogram), as well as ECP navigation both via ‘point & click’ and manual adjustment through incremental tilt and rotations. The GUI has also been written in such a way that features can also be used within text-based scripts, e.g. for repeat matching experiments. To aid in precise analysis of the ECP, it is possible to directly optimize the match of simulation and experiment via the ‘Refine’ button. This refinement algorithm is based on image correlation, by the interior point algorithm as implemented in the MATLAB fmincon function to maximize the normalized cross-correlation between simulated and experimental patterns, similar to EBSD-based pattern matching approaches (Pang et al., 2020). The software also includes a utility to mark a “reference” crystallographic direction and navigate to a new target orientation, with the required stage motion computed and applied in the GUI to align the crystal in the microscope with a particular [uvw] along the optic axis, so that the ECCI/ECP contrast for a new zone axis or a particular band edge can be explored quickly in the SEM. 

-- v1.1 additions --
You can now have access to a refined template matching indexing scheme, based upon the algorithm used in EBSD pattern analysis (Foden et al. 2019, Ultramicroscopy). To use this tool, you should know and set your DD value reasonably well prior to performing the search. This is accessed with the 'index' button. There is also an 'atlas' button which will plot a sterographic projection of the dynamically predicted ECP pattern and the reference sphere, to aid searching for orientations that are beyond the gnomonic projection shown in the simulation tool.

To enable these functions, you need to add the following to the input deck - this is included within the updated decks in this release:
Input_Data.Index=1; %show the Atlas and Index buttons

Full details of this update are found in the updated preprint on the ArXiV.


