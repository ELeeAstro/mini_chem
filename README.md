# Mini-chem

NOTE: This code is in active development and aims to have continual improvements to stability and speed, please report bugs or improvements you find.

## Thermochemical miniature kinetics chemistry module

Mini-chem is a kinetic chemistry network solver primarily for gas giant atmospheric modelling, pared down from the large chemical networks.
This makes use of 'net forward reaction tables', which reduce the number of reactions and species required to be evolved in the ODE solvers significantly.
Mini-chem's NCHO network currently consists of only 13 species with 10 forward reactions (for 20 total including reversed reactions), making it a lightweight and easy to couple network to large scale 3D GCM models, or other models of interest (such as 1D or 2D kinetic modelling efforts).

These are the papers describing mini-chem methods so far: \
Tsai et al. (2022) - A Mini-Chemical Scheme with Net Reactions for 3D GCMs I.: Thermochemical Kinetics - A&A Volume 664, id.A82, 16 pp. - arXiv:2204.04201 \
Lee et al. (2023) - A Mini-Chemical Scheme with Net Reactions for 3D GCMs II.: 3D thermochemical modelling of WASP-39b and HD 189733b - A&A Volume 672, id.A110, 13 pp. - arXiv: arXiv:2302.09525

This repository contains the standalone version of mini-chem, containing the source code and testing modules for users to use for their own purposes.

The src_mini_chem_dlsode is the recommended version (and typically fastest) for most purposes. 
The src_mini_chem_seulex and src_mini_chem_dvode contain seulex and dvode solver standalone variants as an alternative. 
These haven't been checked for threadsafe yet, so only use for non-OpenMP applications, MPI only codes should be ok.

Other directories are currently experimental.

All codes use the [Kahan-Babushka-Neumaier](https://en.wikipedia.org/wiki/Kahan_summation_algorithm) compensated sum algorithm to calculate the sum of the forward and backward rates of each species.
This is generally slower than the (in previous versions) [Pairwise Summation](https://en.wikipedia.org/wiki/Pairwise_summation) method for small timesteps, but is more accurate and actually faster when longer timesteps are considered.
This behavior is also seen in the [FRECKLL](https://ui.adsabs.harvard.edu/abs/2022arXiv220911203A/abstract) kinetic code, which uses a 'distillation' technique to reduce the condition number.
Naive summation typically produces the wrong answer, showing the importance of considering more accurate summation methods.

All source code has been tested with the gcc (gfortran), other compilers are under testing.

## mini_chem.nml

### mini_chem

The namelist that describes the simulation set up:

1. Network - 'HO', 'CHO' or 'NCHO' select network
2. T_in - input temperature [K]
3. P_in - input pressure [Pa]
4. t_step - timestep [s]
5. n_step - number of steps
6. n_sp - number of species (NCHO = 13 including He)
7. data_file - path to _data file
8. sp_file - path to _sp file
9. net_dir - path to net reaction files
10. met - metallicity of net reaction tables

### mini_chem_VMR

1. CE_IC - True or False, interpolate T_in and P_in to CE table in IC file for each mini-chem species.
2. IC_file - Location of initial condition CE file (typically in chem_data/IC).
3. VMR - initial VMR of each species (in the species order of the _sp file).

## Directories

### src_mini_chem_dlsode

Contains a standalone version using the dlsode ODE solver (recommended).

### src_mini_chem_dvode

Contains a standalone version using the dvode ODE solver (alternative).

### src_mini_chem_dlsodes

Contains a standalone version using the dlsodes ODE solver (experimental).

### src_mini_chem_seulex

Contains a standalone version using the seulex ODE solver (alternative).

### src_mini_chem_radau5

Contains a standalone version using the radau5 ODE solver (experimental).

### src_mini_chem_rodas

Contains a standalone version using the rodas ODE solver (experimental).

### src_mini_chem_ros4

Contains a standalone version using the ros4 ODE solver (experimental).

### src_mini_chem_limex

Contains a standalone version using the limex ODE solver (experimental).

### outputs_*

Contains the outputs and plotting scripts for the 0D test code. comp_to_vulcan.py can be used to benchmark against vulcan 0D cases.

### chem_data

Mini-chem formatted data files for each network _sp contains the species data, _data contains the reaction network, 1x, 10x directories etc- contains the netrate tables in mini-chem format for that metallicity and C/O ratio.
IC directory contains chemical equilibrium tables calculated using FastChem for use in the ce interpolation routine, typically used for initial conditions in the GCM. 

### jacobian_conversion

A workspace for converting the Jacobian from VULCAN python to the fortran format - typically is required to be done manually. Feel free to e-mail if new Jacobians are required.

### netrate_contour_plots

Contains python code to produce contour plots of the net forward reaction rate tables.

### netrate_originals

Original VULCAN formatted data for each kinetic network net reaction tables for different metallicity and C/O ratio - typically reformatted using the provided python script.

### vulcan_benchmark_data
 
Contains some benchmarking data from VULCAN.

## Future updates

TODO: mini-photochemistry. \
Improve equilibrium condition detection.

## Compiling

The fortran code can be compiled by altering the makefiles in each src_ directory. 
'make' or 'make debug' (for debug flags) can be used to compile the examples.
The executable is then in the main directory. 
Compiled code can be removed by entering 'make clean'.
You will need to clean and recompile if any changes to the source code are made.


## Coupling guide to GCM and other models

The main.f90 file in each src directory gives an idea of how to interface with mini-chem. Copy the mini_*.f90 and integrator source files into your GCM compilation systems (main.f90 is not needed), then use the mini_ch_ce_interp, mini_ch_i_* and read_react_list modules into your GCM physics routine (see main.f90 for example of this).
Typically, three main subroutines are called for mini-chem operation.

call interp_ce_table(n_sp, T_in, P_in, VMR_IC(:), mu, IC_file)

This takes in a table file and interpolates the T and P to find the VMR and molecular weight [g mol-1] at that temperature and pressure. This is useful for initial conditions, or running a CE scheme or forcing ce if mini-chem fails at certain points.

call read_react_list(data_file, sp_file, net_dir, met)

Reads the reaction list network (see namelist file for examples) - basically set these file paths inside the GCM and call this routine (only) once before the first mini-chem call.

call mini_ch_dlsode(T_in, P_in, t_step, VMR(:), network)

Is the main mini-chem call, takes in a temperature [K], pressure [Pa], time step and current VMR values and network (e.g. 'NCHO'). Call this routine for each GCM cell to perform the kinetic chemistry for that cell. This is usually not called every timestep, but every X times (e.g. 1 hour timesteps or whatever is best for your simulation). The time-step passed to the routine is then the `chemical time-step' for the simulation. For example, if mini-chem is called every second hydrodynamical timestep, then t_step = t_hydro*2.

NOTE: GCM tracers must be in the same order as the _sp file!! i.e. OH must be the first tracer (or passed into mini-chem in VMR index 1 and He the last tracer.)
Helium is assumed to be a passive tracer, so is passed to the routine but is not integrated. 
So you would call mini-chem for a 3D structure like: 

do i \
  do j \
    do k \
      call mini_ch_dlsode(T(i,j,k), P(i,j,k), t_step, q(i,j,k,:), network) \
    end do \
  end do \
end do 

where the last dimension of q is size 13 and is the GCM tracer representing the chemical species VMR in mini-chem.

