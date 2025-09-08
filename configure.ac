# initial information about the project
AC_INIT([IMEX_SfloW2D], [1.0], [mattia.demichielivitturi@ingv.it])

# check if the source folder is available
AC_CONFIG_SRCDIR([src/IMEX_SfloW2D.f90])

# check for fortran compiler
AC_PROG_FC

# automake initialization (mandatory) including a check for automake API version >= 1.9
AM_INIT_AUTOMAKE([1.9 foreign])

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])

AC_ARG_WITH([netcdf],
  AS_HELP_STRING([--with-netcdf=DIR], [location of NetCDF installation]),
  [NETCDF=$withval],
  [NETCDF="/usr/local/netcdf"]
)

AC_SUBST([NETCDF])

AC_PREFIX_DEFAULT([..])

# finally this generates the Makefiles etc. for the build
AC_OUTPUT
