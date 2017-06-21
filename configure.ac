#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(marble, 0.1, hannuksela@phy.cuhk.edu.hk)
AM_INIT_AUTOMAKE([subdir-objects]) # Sub directory cpp and header files
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lcmaes':
AC_CHECK_LIB([cmaes], [main],[HAVE_CMAES=1], [HAVE_CMAES=0])
# FIXME: Replace `main' with a function in `-lcmaes':
AC_CHECK_LIB([openblas], [main],[], [AC_MSG_ERROR(["Could not find openblas"])])
# FIXME: Replace `main' with a function in `-lgsl':
AC_CHECK_LIB([gsl], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
AC_CHECK_LIB([fopenmp], [main])
AC_CHECK_LIB([xerces-c], [main], [HAVE_XERCES_C=1], [HAVE_XERCES_C=0])
AC_CHECK_LIB([gyoto-stdplug], [main],[HAVE_GYOTO_STD=1],[HAVE_GYOTO_STD=0])
AC_CHECK_LIB([gyoto], [main],[HAVE_GYOTO=1], [HAVE_GYOTO=0])
AM_CONDITIONAL([USE_GYOTO], [test "$HAVE_GYOTO" -eq 1])


# Checks for header files.
AC_CHECK_HEADERS([fenv.h float.h limits.h stdlib.h string.h unistd.h])
AC_CHECK_HEADERS([mpi.h],[], [AC_MSG_ERROR([Couldn't find or include mpi.h])])


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 Tests/Makefile
		 PostProcessing/Makefile])
AC_OUTPUT
