# for gfortran Compiler
#===================

F90 = gfortran
F90LINKER =  gfortran

#Debugging and development flags
#FFLAGS	= -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow
#FFLAGS	= -Og -pipe -Wall -Wextra -g -fbacktrace

#Serial flags
FFLAGS  = -O3 -pipe 

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
limex_dense.o \
limex_aux.o \
limex.o \
mini_ch_precision.o \
mini_ch_class.o \
mini_ch_read_reac_list.o \
mini_ch_chem.o \
mini_ch_i_limex.o \
mini_ch_main.o


# executable statement
EXECS  = ../mini_chem_limex

.SUFFIXES : .o .f90 .f .F90

default: mini_chem_mp

mini_chem_mp: $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)
#$(LAPACK) $(BLAS)

clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)

.f90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.F90.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<

.f.o:
	$(F90) $(FFLAGS) $(DEFS) -c $<
