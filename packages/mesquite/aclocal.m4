# Determine which macro to use for the function name.  Possiblities are:
# __PRETTY_FUNCTION__ - g++ - Decorated function name (namespace, class, return type, args, etc.)
# __FUNCDNAME__             - MS Visual Studio decorated function name
# __FUNCTION__ - g++, Visual Studio, others?
# __func__     - C99 standard (not C++), newer g++, (Sun Forte with -features=extensions)
# __function__ - ?
# __FUNC__     - SGI? 
# If __func__ does not work, redefine it to be whichever one does. If
# platform doesn't support any, define it to be an empty string.

AC_DEFUN([MSQ_CPLUSPLUS_FUNC], [
AC_MSG_CHECKING([for function-name preprocessor macro])
AC_LANG_PUSH(C++)
msq_cpp_func=

  # Loop through list of possibilities...
for i in __PRETTY_FUNCTION__ __FUNCDNAME__ __func__ __FUNCTION__ __function__ __FUNC__; do
  if test -z "$msq_cpp_func"; then
    AC_TRY_COMPILE(
      [#include <stdio.h>],
      [(void)printf("%s\n", $i);],
      [msq_cpp_func=$i]
    )
  fi
done

  # If none was found...
if test -z "$msq_cpp_func"; then
  AC_MSG_RESULT(none)
  AC_MSG_WARN([MSQ_FUNCTION will be defined as an empty string.]);
  AC_DEFINE(MSQ_FUNCTION, [""], [Define to c++ preprocessor macro for function name])
else 
  AC_MSG_RESULT($msq_cpp_func)
  AC_DEFINE_UNQUOTED(MSQ_FUNCTION, $msq_cpp_func, [Define to c++ preprocessor macro for function name])
fi

AC_LANG_POP(C++)
]) #MSQ_CPLUSPLUS_FUNC






# Check for existance of new no-file-extension C++ headers.
# Conditinionally sets the following preprocessor macros:
#   MSQ_USE_OLD_STD_HEADERS
#    - The new no-extension versions of the c++ headers do not
#      exist or do not define symbols to be in the "std" namespace.
#      For example: if this is defined #include <map.h> rather than <map>.
#   MSQ_USE_OLD_IO_HEADERS
#    - The new no-extension versions of the c++ IO headers do not
#      exist or do not define symbols to be in the "std" namespace.
#      For example: if this is defined, #include <iostrea.h> rather than <iostream>
#   MSQ_USE_OLD_C_HEADERS 
#    - C++ versions of the standard C headers do not exist or do not place the
#      C symbols in the "std" namespace (e.g. cctype, cstring)
AC_DEFUN([MSQ_CPLUSPLUS_HEADERS], [
AC_LANG_PUSH(C++)

AC_MSG_CHECKING([for C++-standard header files])
AC_TRY_COMPILE([#include <vector>],[std::vector<int> v;],
               [AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no)
                AC_DEFINE(MSQ_USE_OLD_STD_HEADERS,,
                 [Define if old .h headers must be used for C++ standard lib.])])

AC_MSG_CHECKING([for C++-standard I/O header files])
AC_TRY_COMPILE([#include <iostream>],[std::cout << "cout" << std::endl;],
               [AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no)
                AC_DEFINE(MSQ_USE_OLD_IO_HEADERS,,
                 [Define if old .h headers must be used for C++ I/O.])])

AC_MSG_CHECKING([for C++-standard C header files])
AC_TRY_COMPILE([#include <cmath>],[double d = std::sqrt(1.0);],
               [AC_MSG_RESULT(yes)],
               [AC_MSG_RESULT(no)
                AC_DEFINE(MSQ_USE_OLD_C_HEADERS,,
                 [Define if old C .h headers must be used for C++.])])

AC_LANG_POP(C++)
]) # MSQ_CPLUSPLUS_HEADERS






#
# Test if archiver is valid
#

AC_DEFUN([MSQ_TEST_ARCHIVER], [

src_file=conftest.cc
obj_file=conftest.o
lib_file=conftest.a
rm -f $obj_file $lib_file

echo "int main(int argc, char* argv[]) {return 0;}" > $src_file
if $CXX -c $src_file; then
  if ($1 $lib_file $obj_file && ar t $lib_file) > /dev/null 2>&1; then
    m4_ifval([$2], [$2], true)
  m4_ifval([$3], 
  [else
    $3])
  fi
fi

rm -f $src_file $obj_file $lib_file
src_file=
obj_file=
lib_file=
])





#
# Detect C++ archiver
#
AC_DEFUN([MSQ_CPLUSPLUS_ARCHIVER], [
AC_MSG_CHECKING([for C++ archiver])

  # Loop through list of possibilities...
CXX_ARCHIVER=
for i in "$CXX -xar -o" "$CXX -ar -o" "ar cr"; do
  if test -z "$CXX_ARCHIVER"; then
    MSQ_TEST_ARCHIVER( [$i], [CXX_ARCHIVER="$i"] )
  fi
done

if test -z "$CXX_ARCHIVER"; then
  AC_MSG_RESULT([failed])
  AC_MSG_ERROR([No C++ archiver.])
else
  AC_MSG_RESULT([$CXX_ARCHIVER])
  AC_SUBST(CXX_ARCHIVER)
fi

]) # MSQ_CPLUSPLUS_ARCHIVER
