# for gfortran Compiler
#===================


LS_CONFIG  =  -DFULL_ALGEBRA
BLAS=/Users/gl20y334/lapack-3.10.1/librefblas.a
LAPACK=/Users/gl20y334/lapack-3.10.1/liblapack.a

F90 = gfortran
F90LINKER =  gfortran

#Debugging and development flags
#FFLAGS	= -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
#FFLAGS	= -Og -pipe -Wall -Wextra -g -fbacktrace

#Serial flags
FFLAGS  = -O3 -pipe $(LS_CONFIG)

#Parallel flags
#FFLAGS	= -O3 -pipe -fopenmp

# for ifort Compiler
#====================

#F90 = ifort
#F90LINKER = ifort

#FFLAGS   = -O0 -g -traceback -fpp -prec-div -fp-model source -fpe0 -ipo
#FFLAGS    = -O0 -g -traceback -xHost -fpp -fp-model source -qopenmp -ipo
#FFLAGS    = -O3 -xHost -qopenmp -fpp -fp-model source -ipo

#====================

DEFS      =
INCLUDES  =
LFLAGS    = $(FFLAGS)

OBJECTS = \
dc_lapack.o \
LS_Solver.o \
seulex.o \
rodas.o \
radau5.o \
dvode.o \
dlsode.o \
ROS_f90_Integrator.o \
SDIRK_f90_Integrator.o \
mini_ch_precision.o \
mini_ch_class.o \
mini_ch_read_reac_list.o \
mini_ch_chem.o \
mini_ch_i_seulex.o \
mini_ch_i_rodas.o \
mini_ch_i_radau5.o \
mini_ch_i_dvode.o \
mini_ch_i_dlsode.o \
mini_ch_i_Rosenbrock.o \
mini_ch_i_sdirk.o \
mini_ch_main.o \

# lapack.o \
# lapackc.o \
# dlamch.o \


# executable statement
EXECS  = ../mini_chem

.SUFFIXES : .o .f90 .f .F90

default: mini_chem

mini_chem: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS) $(LAPACK) $(BLAS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
