#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [1.0], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE
LDFLAGS="$LDFLAGS -L/opt/local/lib -L/opt/local/lib/gcc49"
CPPFLAGS="$CPPFLAGS -I/opt/local/lib/gcc49/gcc/x86_64-apple-darwin14/4.9.1/include -I/opt/local/include"

LT_INIT
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_PROG_LIBTOOL(libtool)

# Checks for programs.
AC_PROG_CXX([g++-mp-4.9 g++])
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lgomp':
AC_CHECK_LIB([gomp], [main])
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lgslcblas':
AC_CHECK_LIB([gslcblas], [main])

AC_CHECK_LIB([SDL2], [main], [HAS_LIBSDL2=1])
AM_CONDITIONAL([USE_LIBSDL2], [test "$HAS_LIBSDL2" -eq 1])


# Checks for header files.
AC_CHECK_HEADERS([float.h stdint.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT32_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/tumor3d/Makefile
                 src/bin_abc/Makefile
                 src/lunatick/Makefile
                 src/Makefile])
AC_OUTPUT
