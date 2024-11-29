# Compiler settings
#===================

F90 = gfortran
F90LINKER = $(F90)

# Default flags for release build
FFLAGS = -O3 -pipe

# Debugging flags for debugging purposes
DEBUG_FLAGS = -Og -g -pipe -Wall -Wextra -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow

# Include directories and libraries (if any)
INCLUDES =
LFLAGS = $(FFLAGS)
LIBS =

# Source files
OBJECTS = \
    dlsode.o \
    mini_ch_precision.o \
    mini_ch_class.o \
    mini_ch_ce_interp.o \
    mini_ch_read_reac_list.o \
    mini_ch_chem.o \
    mini_ch_i_dlsode.o \
    main.o

# Executable
EXECS = ../mini_chem_dlsode

#====================
# Phony targets
.PHONY: all clean debug release

# Default target (release build by default)
all: $(EXECS)

# Debug build (overrides FFLAGS with debug flags)
debug: FFLAGS = $(DEBUG_FLAGS)
debug: clean all

# Release build (uses default release flags)
release: clean all

# Linking the executable
$(EXECS): $(OBJECTS)
	$(F90LINKER) $(LFLAGS) $(OBJECTS) -o $(EXECS)

# Pattern rule for compiling Fortran 90 files
%.o: %.f90
	$(F90) $(FFLAGS) $(DEFS) -c $<

# Pattern rule for compiling Fortran 77 (.f) files
%.o: %.f
	$(F90) $(FFLAGS) $(DEFS) -c $<

# Clean target (removes object files, module files, and executable)
clean:
	rm -f *.o *.mod *~ *__genmod.f90 $(EXECS)
