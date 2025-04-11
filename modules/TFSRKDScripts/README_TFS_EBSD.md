# TFS_EBSD_Utils

MATLAB Scripts for processing Thermo Fisher RKD/EBSD data (xTalView).

The scripts, especially pattern processing, are based on the work by Jakub Holzer originally developed in pythob.

## MATLAB

The scripts are located in two folders: ``TFSRKDUtils`` (functions) and ``TFSRKDScripts`` (scripts).

### TFSRKDScripts

- ``EBSD_load_data.m``: load orientation data, phase data, and plot maps using [MTEX](https://mtex-toolbox.github.io/index). Need `.cif` files for symmetry.
- ``EBSD_load_processed.m": load processed patterns into a workspace variable (3D array).
  - ``RKD_load_processed.m`` is the RKD equivalent.

- ``bCreateHDF5FromTFSData.m`` (testing): package metadata, orientation data and patterns into a Bruker-like .h5 file.
  - some metadata are loaded using ``ReadTFSJson``
  - Patterns: processed patterns can be loaded (with possible binning) using ``EBSD_load_processed_for_h5.m``. Raw pattern equivalent (testing) is ``EBSD_load_raw_for_h5.m``. RKD equivalents are also available.
  - The equivalent scripts for loading patterns into a MATLAB workspace variable (3D array) is ``EBSD_load_processed.m``. 
- ``SinglePatchEBSP.m``: you can use this script to go through maps in a project with only a single pattern. Useful for parametric studies. Pattern simulation and cross-correlation parts requires [AstroEBSD](https://github.com/ExpMicroMech/AstroEBSD) and appropriate dynamical templates.

### TFSRKDUtils

- ``ReadTFSJson.m``: a wrapper for reading JSON files, since most Thermo Fisher metadata are stored in JSON files.
- ``ReadTFSOrientationData.m``: loading orientation data from ``your_project/results/multiphase.idx`` into a table. 
- ``TileQuadPatterns.m`` (RKD only): tile the four detector quadrants and apply a crop to form a single RKD pattern.