#
# Build all Secondary Stable CASL VRI add-on Trilnos packages with GCC 4.5.1 compiler
#

INCLUDE("${CTEST_SCRIPT_DIRECTORY}/TrilinosCTestDriverCore.pu241.icpc.12.0.4.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/SubmitToCaslDevVera.cmake")
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/casl-vri-packages-coupled.cmake")

SET(COMM_TYPE SERIAL)
SET(BUILD_TYPE DEBUG)
SET(BUILD_DIR_NAME SERIAL_DEBUG_ICPC_CASLDEV)
#SET(CTEST_TEST_TYPE Experimental)
#SET(CTEST_TEST_TIMEOUT 900)
SET(Trilinos_ENABLE_SECONDARY_STABLE_CODE ON)
SET(EXTRA_CONFIGURE_OPTIONS
  # ToDo: Add other extra configure options
  )

TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER()
