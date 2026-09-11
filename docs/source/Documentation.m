%% Documentation
%
%
% *TODO:* ADD Citations properly
% *TODO:* ADD h5 readers and writers that can interact with all the bits
% more nicely
%
%
% In this documentation we TRY to document the many individual skripts we
% developed in analogy to the AstroEBSD repositiory
% Until know there may be functionallity from AstroEBSD used, which is not
% publically available. However, we try to avoid these instances.
% STATUS 2025/08/08 -> Documentation just started. By no means this is
% suppposed to be complete
%
% Main author of decks in astro_skripts: Lukas Berners (many of them are
% taken from the AstroEBSD where main author is T.B. Britton (see Asro EBSD
% repo) Some skripts have contributions from T. Zhang and D. Zeitter.

%% Contents
% The repository contains decks to handle the RTM, PCA and Reprojection
% functionallity from Astro EBSD
%
%


%% RTM
% This contains mostly the function from original Astro Repo.
% We added however, a new optimisation routine which is based again on the
% real space  (using either nelder-mead (see Chris Schuh group publication)
% or interior point minimization (utilised in the AstroECP deck.
% There are multiple ways how we could combine the original Hough based
% indexing wiht the realSpace based RTM and the XCF (Fourier space) based
% RTM
% See the RTM decks for this

%% PCA:
% here we mostly worked on getting the PCA run with EDAX Data and on the
% RWTH HPC computing cluster
% For this we found that it is easier and faster to clone the Astro, MTEX
% and astroskripts repositiory on to the Cluster and run the skript from
% there using a shell script
% *TODO*: Right now the Scripts are cluttered with all the input parameters,
%       try to extract them into individual config files (need to see how this will go)


%% Reprojections
% This Deck Folder contains various reprojections that are based on the
% Gazing at crystall balls paper but borrows some functionallity from the
% Multi Fusion Paper by T. Zhang. 
% This also contains an algorithm that can extract band profiles discretely
% by reprojection of the "crystallographic" band on to the gnomonic EBSD
% patterns and using a gridded interpolant for extraction of the
% intensities


%% Atom plotting
% This Deck contains functionallity to plot atom structures (read in from
% cif)
% here we can Plot the complete structure in 3d, or plot planar sections of
% the crystal structure, where the atoms are colored according to their
% distance to the plane of projection


%% Dynamic Comparisons
% This deck contains functions to compare the dynamically calculated
% templates with each other,
% we utilise this deck to try and find small differences between
% structurally similar arrangements.