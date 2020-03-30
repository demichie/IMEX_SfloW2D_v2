
AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = SW_VAR_DENS_MODEL
bindir=@bindir@

# list of sources for the binary
SW_VAR_DENS_MODEL_SOURCES = \
                parameters_2d.f90 \
                complexify.f90 \
                geometry_2d.f90 \
                constitutive_2d.f90 \
                solver_2d.f90 \
                init_2d.f90 \
                inpout_2d.f90 \
                SW_VAR_DENS_MODEL.f90

FCFLAGS         = -fopenmp -Ofast -funroll-all-loops
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran


