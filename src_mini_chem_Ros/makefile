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
FFLAGS  = -O2 -pipe $(LS_CONFIG)

#Parallel flags

#FFLAGS	= -O0 -g -fbacktrace -pipe -fopenmp
#FFLAGS	= -O3 -pipe -fopenmp


# for ifort Compiler
#====================

#F90 = ifort
#F90LINKER = ifort

#FFLAGS   = -O0 -g -traceback -fpp -prec-div -fp-model source -fpe0 -ipo
#FFLAGS    = -O0 -g -traceback -xHost -fpp -fp-model source -qopenmp -ipo
#FFLAGS    = -O3 -qopenmp -fpp -fp-model source -ipo

#====================

DEFS      =
INCLUDES  =
LFLAGS    = $(FFLAGS) 

OBJECTS = \
LS_Solver.o \
ROS_f90_Integrator.o \
mini_ch_precision.o \
mini_ch_class.o \
mini_ch_read_reac_list.o \
mini_ch_chem.o \
mini_ch_i_Ros.o \
mini_ch_main.o


# executable statement
EXECS  = ../mini_chem_Ros

.SUFFIXES : .o .f90 .f .F90

default: mini_chem_Ros

mini_chem_Ros: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS) $(LAPACK) $(BLAS)
#$(LAPACK) $(BLAS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
