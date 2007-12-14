
CONFIGURE_FILE(${RBGen_SOURCE_DIR}/cmake/RBGen_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/RBGen_config.h)

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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${Teuchos_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${Teuchos_BINARY_DIR}/src)

ADD_LIBRARY(rbgen ${HEADERS} ${SOURCES})

TARGET_LINK_LIBRARIES(rbgen teuchos)

