#------------------------------------------------
# Experimental LINUX tests
# Debug with Coverage and MemoryCheck
#------------------------------------------------
SET (CTEST_SOURCE_NAME Trilinos)
SET (TEST_TYPE nightly)
SET (BUILD_TYPE debug)

SET (CTEST_DASHBOARD_ROOT /var/dashboards)
SET (CTEST_CMAKE_COMMAND "\"${CMAKE_EXECUTABLE_NAME}\"")

# Options for Nightly builds
SET (CTEST_BACKUP_AND_RESTORE TRUE)
SET (CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)
SET (CTEST_CVS_CHECKOUT
  "cvs -Q -d :ext:$ENV{USER}@software.sandia.gov:/space/CVS co -d \"${CTEST_SOURCE_NAME}\" ${CTEST_SOURCE_NAME}"
)

SET (CTEST_BINARY_NAME Trilinos-${TEST_TYPE}-${BUILD_TYPE})
SET (CTEST_SOURCE_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_SOURCE_NAME}")
SET (CTEST_BINARY_DIRECTORY "${CTEST_DASHBOARD_ROOT}/${CTEST_BINARY_NAME}")

SET (CTEST_COMMAND 
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalStart"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalUpdate"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalConfigure"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalBuild"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalSubmit"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalTest"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalCoverage"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalSubmit"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalMemCheck"
  "\"${CTEST_EXECUTABLE_NAME}\" -D ExperimentalSubmit"
)

SET (CTEST_INITIAL_CACHE "

TRILINOS_ENABLE_ANASAZI:BOOL=ON
TRILINOS_ENABLE_EPETRA:BOOL=ON
TRILINOS_ENABLE_RBGEN:BOOL=ON
TRILINOS_ENABLE_TEUCHOS:BOOL=ON

TRILINOS_ENABLE_MPI:BOOL=ON

EPETRA_ENABLE_FORTRAN:BOOL=ON
RBGEN_ENABLE_EPETRA:BOOL=ON
TEUCHOS_ENABLE_COMPLEX:BOOL=ON
TEUCHOS_ENABLE_EXTENDED:BOOL=ON

TEUCHOS_ENABLE_TESTS:BOOL=ON
ANASAZI_ENABLE_TESTS:BOOL=ON
EPETRA_ENABLE_TESTS:BOOL=ON

TEUCHOS_ENABLE_EXAMPLES:BOOL=OFF
EPETRA_ENABLE_EXAMPLES:BOOL=OFF
ANASAZI_ENABLE_EXAMPLES:BOOL=OFF

BUILDNAME:STRING=$ENV{HOSTTYPE}-${TEST_TYPE}-${BUILD_TYPE}

CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}

CMAKE_CXX_FLAGS:STRING=-g -O0 -Wall -W -Wshadow -Wunused-variable -Wunused-function -Wno-system-headers -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage

CMAKE_C_FLAGS:STRING=-g -O0 -Wall -W -fprofile-arcs -ftest-coverage

CMAKE_EXE_LINKER_FLAGS:STRING=-fprofile-arcs -ftest-coverage

MAKECOMMAND:STRING=gmake -j 4
")


