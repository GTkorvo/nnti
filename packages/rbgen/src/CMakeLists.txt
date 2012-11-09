
INCLUDE(TribitsLibraryMacros)
INCLUDE(TribitsAddExecutable)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
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

SET(SOURCES ${SOURCES}
  RBGen_Params.cpp
  )

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
IF(${PACKAGE_NAME}_ENABLE_Epetra)

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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF(${PACKAGE_NAME}_ENABLE_EpetraExt)

  SET(HEADERS
    ${HEADERS}
    RBGen_EpetraCrsMatrixFileIOHandler.h
    RBGen_MatrixMarketFileIOHandler.h
    )

  SET(SOURCES
    ${SOURCES}
    RBGen_EpetraCrsMatrixFileIOHandler.cpp
    RBGen_MatrixMarketFileIOHandler.cpp
  )

ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Anasazi)
IF(${PACKAGE_NAME}_ENABLE_Anasazi AND ${PACKAGE_NAME}_ENABLE_Epetra)

  SET(HEADERS
    ${HEADERS}
    RBGen_AnasaziPOD.h
    RBGen_IncSVDPOD.h
    RBGen_ISVDMultiCD.h
    RBGen_ISVD_MultiCDUDV.h
    RBGen_ISVDMultiSDA.h
    RBGen_ISVD_MultiSDAUDV.h
    RBGen_ISVDMultiSDB.h
    RBGen_ISVD_MultiSDBUDV.h
    RBGen_ISVDSingle.h
    RBGen_ISVD_SingleUDV.h
    RBGen_ISVDUDV.h
    RBGen_StSVD_RTR.h
    )

  SET(SOURCES
    ${SOURCES}
    RBGen_AnasaziPOD.cpp
    RBGen_IncSVDPOD.cpp
    RBGen_ISVDMultiCD.cpp
    RBGen_ISVD_MultiCDUDV.cpp
    RBGen_ISVDMultiSDA.cpp
    RBGen_ISVD_MultiSDAUDV.cpp
    RBGen_ISVDMultiSDB.cpp
    RBGen_ISVD_MultiSDBUDV.cpp
    RBGen_ISVDSingle.cpp
    RBGen_ISVD_SingleUDV.cpp
    RBGen_ISVDUDV.cpp
    RBGen_StSVD_RTR.cpp
  )

ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  rbgen
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Define the targets for a driver
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Epetra)
IF(${PACKAGE_NAME}_ENABLE_Epetra)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF(${PACKAGE_NAME}_ENABLE_EpetraExt)

  SET(EPETRADRIVER_SOURCES
    ${SOURCES}
    RBGenDriver_EpetraMV.cpp
    )

  TRIBITS_ADD_EXECUTABLE(
    RBGenDriver_EpetraMV
    NOEXEPREFIX
    SOURCES ${EPETRADRIVER_SOURCES}
    COMM serial mpi
  )

ENDIF()

ENDIF()

