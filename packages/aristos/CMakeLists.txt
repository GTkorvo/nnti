

#
# A) Define the package
#

TRIBITS_PACKAGE(Aristos)

#
# B) Set up package-specific options
#

IF(Amesos_ENABLE_UMFPACK)
  SET(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK_DEFAULT OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_AMESOS_UMFPACK
  HAVE_ARISTOS_AMESOS_UMFPACK
  "Enable Amesos with UmfPack functionality."
  ${Aristos_ENABLE_AMESOS_UMFPACK_DEFAULT} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_ARISTOS_DEBUG
  "Enable a host of runtime debug checking."
  ${Trilinos_ENABLE_DEBUG} )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
