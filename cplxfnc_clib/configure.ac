AC_PREREQ([2.69])
AC_INIT([cplxfnc],[0.1],[https://github.com/cimatosa/cplxfnc/issues])
AC_CONFIG_SRCDIR([cplxfnc.hpp])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
#AC_PROG_INSTALL
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [exp],,AC_MSG_ERROR([Unable to find math library (libm)]))
AC_CHECK_LIB([flint], [fmpz_init],,AC_MSG_ERROR([Unable to find flint library (libflint)]))
#AC_SEARCH_LIBS([acb_hurwitz_zeta],[arb flint-arb],,AC_MSG_ERROR([Unable to find arb library (libarb or libflint-arb)]))

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.

AC_CONFIG_FILES([makefile])
AC_OUTPUT
