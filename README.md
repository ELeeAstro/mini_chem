# Mini-chem

## Thermochemical miniature kinetics chemistry module

Mini-chem is a kinetic chemistry network solver primarily for gas giant atmospheric modelling, pared down from the large chemical networks.
This makes use of 'net forward reaction tables', which reduce the number of reactions and species required to be evolved in the ODE solvers significantly.
Mini-chem's NCHO network currently consists of only 12 species with 10 reactions, making it a lightweight and easy to couple network to large scale 3D GCM models, or other models of interest (such as 1D or 2D kinetic modelling efforts).

These are the papers describing mini-chem methods so far: \
Tsai et al. (2022) - A Mini-Chemical Scheme with Net Reactions for 3D GCMs I.: Thermochemical Kinetics - A&A Volume 664, id.A82, 16 pp. - arXiv:2204.04201 \
Lee et al. (2023) - A Mini-Chemical Scheme with Net Reactions for 3D GCMs II.: 3D thermochemical modelling of WASP-39b and HD 189733b - A&A  Volume 672, id.A110, 13 pp. - arXiv: arXiv:2302.09525

This repository contains the standalone version of mini-chem, containing the source code and testing modules for users to use for their own purposes.

The src_mini_chem_seulex is the reccomended version for most purposes. 
The src_mini_chem_dvode, src_mini_chem_dlsode and src_mini_chem_dlsodes contain dvode, dlsode and dlsodes solver standalone variants as an alternative. 
These haven't been checked for threadsafe yet, so only use for non-OpenMP applications.

Other directories are currently experimental

NOTE: Integration of the network may take a long time for high pressures (> 10-20 bar) and high temperatures (> 1000 K), we suggest applying a CE scheme for those deep regions, which should be at CE anyway.

This code is in active development and aims to have continual improvements to stability and speed, please report bugs or improvements you find.

NOTE: The rodas solver may not work with intel compilers, it is unclear why. 

## mini_chem.nml

### mini_chem

The namelist that describes the simulation set up:

1. Network - 'HO', 'CHO' or 'NCHO' select network
2. T_in - input temperature [K]
3. P_in - input pressure [Pa]
4. t_step - timestep [s]
4. n_step - number of steps
5. data_file - path to _data file
6. sp_file - path to _sp file
7. net_dir - path to net reaction files
8. met - metallicity of net reaction tables

### mini_chem_VMR

1. VMR - initial VMR of each species (in the species order of the _sp file)

## Directories

### src_mini_chem_seulex

The standalone seulex version of mini_chem, reccomended for most applications. 

### src_mini_chem_mp

Experimental

### src_mini_chem_dvode

Contains a standalone version using the dvode ODE solver (alternative)

### src_mini_chem_dlsode

Contains a standalone version using the dlsode ODE solver (alternative)

### src_mini_chem_dlsodes

Contains a standalone version using the dlsodes ODE solver (alternative)

### chem_data

Mini-chem formatted data files for each network _sp contains the species data, _data contains the reaction network, 1x, 10x directories etc- contains the netrate tables in mini-chem format for that metallicity and C/O ratio

### jacobian_conversion

A workspace for converting the Jacobian from VULCAN to a FORTRAN format - typically is required to be done manually. Feel free to e-mail if new Jacobians are required

### netrate_contour_plots

Contains python code to produce contour plots of the net forward reaction rate tables

### netrate_originals

Original VULCAN formatted data for each kinetic network net reaction tables for different metallicity and C/O ratio - typically reformatted using the provided python script

### outputs, outputs_mp, outputs_dvode etc

Where the output from running the code is produced and some python plotting routines.

### vulcan_benchmark_data
 
Contains some benchmarking data from VULCAN.

## Future updates

TODO: Photochemistry and haze \
Improve efficency of pairwise addition implimentation \
Improve equilibrium condition detection

## Compiling

The fortran code can be compiled by altering the makefiles in each src_ directory.
Compiled code can be removed by entering 'make clean'.
