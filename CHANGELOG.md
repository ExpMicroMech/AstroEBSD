# Change Log

Here we list major updates to the repository since 25/09/2019. Small fixes are listed under the corresponding major versions.

Prof. Britton (TBB) is involved in most of these updates, and additional contributers are also listed. 

## 10/07/2025 - AstroECP

See `modules/AstroECP` for details.

- Added AstroECP as a pattern indexing and visualization tool and an experimental control assistant that facilitates ECCI by precisely dictating about the requisite stage tilts/rotations of SEM to navigate along the crystal.
- A digital twin of our SEM within a graphical user interface is included to help navigate the electron channeling patterns and quickly access specific ECCI conditions for subsequent analysis.

- Edited the main `readme` file, `credit`, and added `CHANGELOG`.

Principal contributions from M. Haroon Qaiser

## 08/07/2025 - Thermo Fisher data loader (v1.4.0)

See `modules/TFSRKDScripts` for details

- Added code for loading data from Thermo Fisher Scientific systems, including TruePix (EBSD) and RKD data. This include:
  - Reading and plotting orientation and phase data
  - Reading and plotting other map data, including virtual diode image data
  - Reading raw and background corrected patterns
  - Packaging data into Bruker-like .h5 files 

With help from Tianbi Zhang, Lukas Berners

## 18/02/2025 - Kikuchi geometry comparison (v1.3.0)

General: 

- `phases/masterpattern` renamed to `phases/dynamical_templates`
- Added support of dynamical templates from EMsoft (.h5), AztecCrystal Mapsweeper and the Winkelmann tool ("BWKD"). A few examples of Si are provided.
- A tutorial document on how to generate dynamical templates and use them in AstroEBSD is added to `tutorial_docs`.

Kikuchi geometry comparison (see `modules/ded_geometry_comparison`)

- Improved code for Spherical reprojection in Cartesian coordinates
- Codes for spherical reprojection in spherical coordinates and band profile analysis
- Codes for equivalent pattern quality comparison
- Miscellaneous updates of phase files and .cif files.

With help from Tianbi Zhang, Lukas Berners

## 02/10/2024 - new h5oina loader

With help from Dr. Ruth Birch

## 23/11/2023 - direct electron detector (MiniPIX TPX3)-based modular systems (v1.2.0)

See `modules/ded-tkd`, `modules/ded-ebsd` for details.

- Codes for exposure fusion of on-axis TKD patterns
- Geometry calibration and reindexing of static DED-based EBSD system
- Spherical reprojection in Cartesian coordinates

With help from Tianbi Zhang


## 12/02/2020 - Tom McAuliffe

- New phase data structures: phase files, cifnames and masterpatterns (.bin files) specified. 
- Phase_Builder_RTM function updated accordingly to handle these.

- Integration of the RTM and PCA modules.
  - RTM example deck for indexing of a full map with a single candidate phase.
  - PCA example decks for clustering of combined EBSD/EDS datasets, and EBSD-based phase ID and indexing.
  - PCA postanalysis decks for chemical analysis (as function of phase) and RC-EBSP comparison to template simulations.
  - Available with (adjustable) spatial weighting kernel.

- .gitignore for phase files / folders fixed


## 25/09/2019

- cyan/yellow colouring fix

- new background correction options

- mac/windows loader fixed

With help from Alex Foden & Tom McAuliffe

