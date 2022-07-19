# Mini-chem

## Thermochemical miniature kinetics chemistry module

Mini-chem is a kinetic chemistry network solver primarily for gas giant atmospheric modelling, pared down from the large chemical networks.
This makes use of 'net forward reaction tables', which reduce the number of reactions and species required to be evolved in the ODE solvers significantly.
Mini-chem's NCHO network current consists of only 12 species with 10 reactions, making it a lightweight and easy to couple network to large scale 3D GCM models, or other models of interest (such as 1D or 2D kinetic modelling efforts).

These are the papers describing our methods so far: \
Tsai et al. (2022) - A Mini-Chemical Scheme with Net Reactions for 3D GCMs I.: Thermochemical Kinetics - A&A (in review) - arXiv:2204.04201 \
Lee et al. (202X) - A Mini-Chemical Scheme with Net Reactions for 3D GCMs IB.: 3D thermochemical modelling of WASP-39b and HD 189733b - A&A (submitted)

This repository contains the standalone version of mini-chem, containing the source code and testing modules for users to use for their own purposes.

The branch 'multi_threading' contains the threadsafe code, suitable for GCM modelling. We recommend using the source files from this branch as the production methods inside the GCM model.

## mini_chem.nml

### mini_chem

The namelist that describes the simulation set up:

1. Network - 'HO', 'CHO' or 'NCHO' select network
2. T_in - input temperature [K]
3. P_in - input pressure [Pa]
4. t_step - timestep
4. n_step - number of steps
5. data_file - path to _data file
6. sp_file - path to _sp file

### mini_chem_VMR

1. VMR - initial VMR of each species (in the species order of the _sp file)


## Directories

### src_mini_chem

Contains the fortran source files and makefile for compiling this standalone version

### chem_data

Mini-chem formatted data files for each network _sp contains the species data, _data contains the reaction network, solar- contains the netrate tables in mini-chem format

### jacobian_conversion

A workspace for converting the Jacobian from VULCAN to FORTRAN - typically is required to be done manually

### netrate_contour_plots

Contains python code to produce contour plots of the net forward reaction rate tables

### netrate_originals

Original VULCAN formatted data for each kinetic network net reaction tables - typically reformatted using the python script

## outputs

Where the output from running the code is produced, also contains some benchmarking data from VULCAN

### Future updates

We will be adding photochemical effects in the next projects. \
Add Bezier interpolation to net rate table interpolation for multi_threading. \
