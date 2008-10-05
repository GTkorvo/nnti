
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PROJECT_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS
  RBGen_ConfigDefs.h 
  RBGen_FileIOFactory.hpp 
  RBGen_FileIOHandler.hpp 
  RBGen_Filter.hpp 
  RBGen_MethodFactory.hpp 
  RBGen_Method.hpp 
  RBGen_NoPreprocessor.hpp 
  RBGen_Params.h 
  RBGen_PODMethod.hpp 
  RBGen_PreprocessorFactory.hpp 
  RBGen_Preprocessor.hpp
  )

SET(SOURCES
  RBGen_Params.cpp
  )

IF(${PROJECT_NAME}_ENABLE_Epetra)

  SET(HEADERS
    ${HEADERS}
    RBGen_BurkardtFileIOHandler.h 
    RBGen_EpetraMVFileIOFactory.h 
    RBGen_EpetraMVMethodFactory.h 
    RBGen_EpetraMVPreprocessorFactory.h 
    RBGen_LapackPOD.h 
    RBGen_MSPreprocessor.h 
    RBGen_Utils.h
    )

  SET(SOURCES
    ${SOURCES}
    RBGen_BurkardtFileIOHandler.cpp 
    RBGen_EpetraMVFileIOFactory.cpp 
    RBGen_EpetraMVMethodFactory.cpp 
    RBGen_EpetraMVPreprocessorFactory.cpp 
    RBGen_LapackPOD.cpp 
    RBGen_MSPreprocessor.cpp 
    RBGen_Utils.cpp
    )

ENDIF()


#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  rbgen
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
