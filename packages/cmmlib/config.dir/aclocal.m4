dnl 
dnl modified by Mark Sears
dnl I deleted fortran stuff, compacted all into aclocal.m4
dnl

dnl @synopsis MPS_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl Example use
dnl
dnl MPS_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1],
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])

dnl @synopsis MPS_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl
dnl Prepends the specified name(s) to the list of libraries to link
dnl with.
dnl
dnl Example use
dnl
dnl MPS_ARG_WITH_LIBS
dnl
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([MPS_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs],
[List additional libraries here.  For example, --with-libs=-lsuperlu or --with-libs=/path/libsuperlu.a]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)

dnl @synopsis MPS_ARG_WITH_LIBDIRS
dnl
dnl Test for --with-libdirs="-Llibdir1 -Llibdir2".  if defined,
dnl prepend "-Llibdir1 -Llibdir2" to LDFLAGS
dnl
dnl Use this macro to facilitate addition of directories to library search path.
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_WITH_LIBDIRS],
[
AC_MSG_CHECKING([whether additional library search paths defined])
AC_ARG_WITH(libdirs,
AC_HELP_STRING([--with-libdirs],
[OBSOLETE use --with-ldflags instead.  (ex. --with-ldflags="-L<DIR> -L<DIR2>")]),
[
LDFLAGS="${withval} ${LDFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis MPS_ARG_WITH_INCDIRS
dnl
dnl Test for --with-incdirs="-Iincdir1 -Iincdir2".  if defined, prepend
dnl "-Iincdir1 -Iincdir2" to CPPFLAGS
dnl
dnl Use this macro to facilitate addition of directories to include file search path.
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_WITH_INCDIRS],
[
AC_MSG_CHECKING([whether additional include search paths defined])
AC_ARG_WITH(incdirs,
AC_HELP_STRING([--with-incdirs],
[additional directories containing include files: will prepend to search here for includes, use -Idir format]),
[
CPPFLAGS="${withval} ${CPPFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])

dnl @synopsis MPS_ARG_WITH_AR
dnl
dnl Test for --with-ar="ar_program ar_flags".
dnl Default is "ar cru"
dnl
dnl Generates an Automake conditional USE_ALTERNATE_AR that can be tested.
dnl Generates the user-specified archiver command in @ALTERNATE_AR@.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_WITH_AR],
[
AC_ARG_WITH(ar,
AC_HELP_STRING([--with-ar], [override archiver command (default is "ar cru")]),
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

if test -n "${SPECIAL_AR}" && test "X${USE_ALTERNATE_AR}" != "Xyes";
then
  USE_ALTERNATE_AR=yes
  ALTERNATE_AR="${SPECIAL_AR}"
fi

AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
#   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
#   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)
])

dnl @synopsis MPS_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --enable-mpi       - Turns MPI compiling mode on
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-compilers - Turns on MPI compiling mode and sets the MPI C++
dnl                       compiler = mpicxx or mpiCC (if mpicxx not available),
dnl                       the MPI C compiler = mpicc and 
dnl                       the MPI Fortran compiler = mpif77
dnl --with-mpi-incdir - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl If any of these options are set, HAVE_MPI will be defined for both
dnl Autoconf and Automake, and HAVE_MPI will be defined in the
dnl generated config.h file
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_CONFIG_MPI],
[

AC_ARG_ENABLE(mpi,
[AC_HELP_STRING([--enable-mpi],[MPI support])],
[HAVE_PKG_MPI=$enableval],
[HAVE_PKG_MPI=no]
)

AC_ARG_WITH(mpi-compilers,
[AC_HELP_STRING([--with-mpi-compilers=PATH],
[use MPI compilers mpicc, mpif77, and mpicxx (or mpiCC) in the specified path or in the default path if no path is specified. Enables MPI])],
[
  if test X${withval} != Xno; then
    HAVE_PKG_MPI=yes
    if test X${withval} = Xyes; then
      # Check for mpicxx, if it does not exist, use mpiCC instead.
      AC_CHECK_PROG(MPI_CXX, mpicxx, mpicxx, mpiCC)
      MPI_CC=mpicc
#      MPI_F77=mpif77
    else
      MPI_TEMP_CXX=${withval}/mpicxx
      if test -f ${MPI_TEMP_CXX}; then
        MPI_CXX=${MPI_TEMP_CXX}
      else
        MPI_CXX=${withval}/mpiCC
      fi
      MPI_CC=${withval}/mpicc
#      MPI_F77=${withval}/mpif77
    fi
  fi
]
)

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory (enables MPI)])],
[
  HAVE_PKG_MPI=yes
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)

#AC_ARG_WITH(mpi-include,
#[AC_HELP_STRING([--with-mpi-include],[Obsolete.  Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.])],
#[AC_MSG_ERROR([--with-mpi-include is an obsolte option.   Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.  For example '--with-mpi-incdir=/usr/lam_path/include'.])]
#)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(HAVE_MPI,,[define if we want to use MPI])
fi

dnl Define Automake version of HAVE_MPI if appropriate

#AM_CONDITIONAL(HAVE_MPI, [test "X${HAVE_PKG_MPI}" = "Xyes"])

dnl
dnl --------------------------------------------------------------------
dnl Check for MPI compilers (must be done *before* AC_PROG_CXX,
dnl AC_PROG_CC and AC_PROG_F77)
dnl 
dnl --------------------------------------------------------------------

if test -n "${MPI_CXX}"; then
  if test -f ${MPI_CXX}; then
    MPI_CXX_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  fi

  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
    CXX=${MPI_CXX}
  else
    echo "-----"
    echo "Cannot find MPI C++ compiler ${MPI_CXX}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C++ compiler using CXX=<compiler>"
    echo "Do not use --with-mpi-compilers if using CXX=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

if test -n "${MPI_CC}"; then
  if test -f ${MPI_CC}; then
    MPI_CC_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
  fi

  if test "X${MPI_CC_EXISTS}" = "Xyes"; then
    echo "Using MPI compiler ${MPI_CC}"
    CC=${MPI_CC}
  else
    echo "-----"
    echo "Cannot find MPI C compiler ${MPI_CC}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C compiler using CC=<compiler>"
    echo "Do not use --with-mpi-compilers if using CC=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C compiler (${MPI_CC}) not found.])
  fi
fi

if test -n "${MPI_F77}"; then
  if test -f ${MPI_F77}; then
    MPI_F77_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
  fi

  if test "X${MPI_F77_EXISTS}" = "Xyes"; then
    F77=${MPI_F77}
  else
    echo "-----"
    echo "Cannot find MPI Fortran compiler ${MPI_F77}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a Fortran 77 compiler using F77=<compiler>"
    echo "Do not use --with-mpi-compilers if using F77=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI Fortran 77 compiler (${MPI_F77}) not found.])
  fi
fi
])

dnl @synopsis MPS_ARG_CHECK_MPI
dnl
dnl Check to make sure any definitions set in MPS_ARG_CONFIG_MPI
dnl are valid, set the MPI flags.  Test MPI compile using C++ compiler.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_CHECK_MPI],
[

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
  fi

  AC_LANG_CPLUSPLUS
  AC_MSG_CHECKING(for mpi.h)
  AC_TRY_CPP([#include "mpi.h"],
    [AC_MSG_RESULT(yes)],
    [
     AC_MSG_RESULT(no)
     echo "-----"
     echo "Cannot link simple MPI program."
     echo "Try --with-mpi-compilers to specify MPI compilers."
     echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
     echo "to specify all the specific MPI compile options."
     echo "-----"
     AC_MSG_ERROR(MPI cannot link)
    ])


  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS="${MPI_LIBS} ${LIBS}"
  fi

#   AC_LANG_CPLUSPLUS
#   AC_MSG_CHECKING(whether MPI will link using C++ compiler)
#   AC_TRY_LINK([#include <mpi.h>],
#   [int c; char** v; MPI_Init(&c,&v);],
#   [AC_MSG_RESULT(yes)],
#   [AC_MSG_RESULT(no)
#    echo "-----"
#    echo "Cannot link simple MPI program."
#    echo "Try --with-mpi-cxx to specify MPI C++ compile script."
#    echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
#    echo "to specify all the specific MPI compile options."
#    echo "-----"
#    AC_MSG_ERROR(MPI cannot link)]
#   )

fi
])


dnl @synopsis MPS_ARG_WITH_BLASLIB
dnl
dnl Test for --with-blaslib="name".
dnl MPS add test for --with-blasinclude="stuff"
dnl
dnl Prepends the specified name to the list of files to check for BLAS
dnl routines.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([MPS_ARG_WITH_BLASLIB],
[
 AC_ARG_WITH(blasinclude,
    AC_HELP_STRING([--with-blasinclude], [definition of BLAS_INCLUDE]),
    [
     BLAS_INCLUDE=${withval}
    ]
  )

 AC_ARG_WITH(blaslib,
    AC_HELP_STRING([--with-blaslib], [name of library containing BLAS: will search lib directories for -lname]),
    [
     USE_BLASLIB=yes
     NEWBLASLIB=${withval}
    ]
  )

 BLASLIBS="cxml blas complib.sgimath"

 if test "X${USE_BLASLIB}" = "Xyes"; then
   BLASLIBS="${NEWBLASLIB} ${BLASLIBS}"
   echo "BLASLIBS ${BLASLIBS}"
 fi
])


dnl @synopsis MPS_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by MPS_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id$
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([MPS_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl MPS_BLAS
