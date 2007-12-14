
CONFIGURE_FILE(${RBGen_SOURCE_DIR}/cmake/RBGen_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/RBGen_config.h)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/RBGen_config.h DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR})

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

IF(RBGEN_ENABLE_EPETRA)
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
ENDIF(RBGEN_ENABLE_EPETRA)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)

IF(RBGEN_ENABLE_EPETRA)
  INCLUDE_DIRECTORIES(${Epetra_SOURCE_DIR}/src)
  INCLUDE_DIRECTORIES(${Epetra_BINARY_DIR}/src)
ENDIF(RBGEN_ENABLE_EPETRA)

ADD_LIBRARY(rbgen ${HEADERS} ${SOURCES})

TARGET_LINK_LIBRARIES(rbgen teuchos)

IF(RBGEN_ENABLE_EPETRA)
  TARGET_LINK_LIBRARIES(rbgen epetra)
ENDIF(RBGEN_ENABLE_EPETRA)

INSTALL(TARGETS rbgen
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
	)

INSTALL(FILES ${HEADERS} DESTINATION ${TRILINOS_INSTALL_INCLUDE_DIR})

