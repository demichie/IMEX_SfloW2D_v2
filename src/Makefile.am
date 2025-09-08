
AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = IMEX_SfloW2D
bindir=@bindir@

# list of sources for the binary
IMEX_SfloW2D_SOURCES = \
                parameters_2d.f90 \
                complexify.f90 \
                geometry_2d.f90 \
                constitutive_2d.f90 \
                solver_2d.f90 \
                init_2d.f90 \
                stochastic_module.f90 \
                inpout_2d.f90 \
                IMEX_SfloW2D.f90

FCFLAGS         = -fopenmp -Ofast -funroll-all-loops -I@NETCDF@/include -lnetcdff
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack -L@NETCDF@/lib -lnetcdff 

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran
NETCDF = @NETCDF@


