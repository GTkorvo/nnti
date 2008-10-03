#
# A) Define the package
#

INCLUDE(TrilinosPackageHelpers)
TRILINOS_PACKAGE(RBGen)

SET(${PACKAGE_NAME}_VERSION "0.1d")

#
# B) Set up package-specific options
#

OPTION(${PACKAGE_NAME}_ENABLE_EPETRA
  "Enable ${PACKAGE_NAME} features that use epetra."
  ON
  )

# 2008/10/03: rabartl: ToDo: Above, the option RBGen_ENABLE_EPETRA should get
# added automatically given the package dependency
# RBGen_LIB_OPTIONAL_DEP_PACKAGES

#
# C) Add the libraries and examples
#

ADD_SUBDIRECTORY(src)

#IF(${PACKAGE_NAME}_ENABLE_TESTS)
#  ADD_SUBDIRECTORY(test)
#ENDIF()
#
#IF(${PACKAGE_NAME}_ENABLE_EXAMPLES)
#  ADD_SUBDIRECTORY(example)
#ENDIF()

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
