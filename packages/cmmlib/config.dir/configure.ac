
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# note that macros beginning with TAC_ are Trilinos macros

# @HEADER
#  See the file "COPYRIGHT" in the CMMLIB distribution.
# @HEADER

#  #np# - package name, version number, and e-mail address below
AC_INIT(cmmlib, 1.0, mpsears@sandia.gov)

echo "----------------------------------------"
echo "Running CMMLIB Configure Script"
echo "----------------------------------------"

# specify required version of autoconf
AC_PREREQ(2.59)

# AC_CONFIG_AUX_DIR(config.dir)

# protection
AC_CONFIG_SRCDIR(./cmmlib.mk.in)

# not sure what this does
AC_CANONICAL_TARGET

# MPS_ macros are put into the aclocal.m4 file.
# I rebuilt these macros because some of them depended on AM_...
# CONFIG_MPI must be checked before AC_PROG_CC etc
MPS_ARG_CONFIG_MPI

# Checks for programs (compilers, etc)
AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)
AC_PROG_RANLIB

# ------------------------------------------------------------------------
# Checks for extra flags
# #np# - You can pull flag checks for languages not used in package.  Please
#         note that if you use F77_FUNC to determine Fortran name mangling,
#         you should not remove the Fortran compiler check or the check for
#         Fortran flags.  Doing so will prevent the detection of the proper
#         name mangling in some cases.
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

MPS_ARG_WITH_FLAGS(ccflags, CCFLAGS)
MPS_ARG_WITH_FLAGS(cflags, CFLAGS)
MPS_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
MPS_ARG_WITH_LIBS
MPS_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------
MPS_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------
echo "CHECK MPI!"
MPS_ARG_CHECK_MPI

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
#TAC_ARG_WITH_LIBDIRS

# enable use of --with-incdirs="-Lincdir1 -Lincdir2 ..." to prepend to CPPFLAGS
#TAC_ARG_WITH_INCDIRS

# check standard C math libs
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR([Cannot find math library],1))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR([Cannot find math library],1))

# Check for BLAS
echo "CHECK BLAS!"
MPS_ARG_WITH_BLASLIB

echo "USE_BLAS ${USE_BLAS}"
echo "NEWBLASLIB ${NEWBLASLIB}"

MPS_BLAS([], AC_MSG_ERROR([Cannot find blas library]))

#AC_MSG_ERROR([Cannot find blas library, specify a path using --with-ldflags=-I<DIR> 
#(ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB 
#(for example  --with-blas=/usr/path/lib/libcxml.a)],1))

AC_SUBST(BLAS_INCLUDE)
AC_SUBST(sgemm)

#
# Specify files to generate (Makefile.in -> Makefile, etc.)
#
AC_CONFIG_FILES([
                package.mk
                cmmlib.mk
                config.h
                ])

AC_OUTPUT()

echo "----------------------------------------"
echo "Finished Running CMMLIB Configure Script"
echo "----------------------------------------"
