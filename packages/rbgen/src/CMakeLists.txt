# $Header$ 

# Setup direct dependencies
SET(DIRECT_DEPS Teuchos)
IF(${PROJECT_NAME}_ENABLE_EPETRA)
  LIST(APPEND DIRECT_DEPS Epetra)
ENDIF(${PROJECT_NAME}_ENABLE_EPETRA)

# Create <PACKAGE>_config.h
CONFIGURE_FILE(
  ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}_config.h.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
)
INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
  DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Setup directories based on current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Setup directories based on dependency directories
FOREACH(DIRECT_DEP ${DIRECT_DEPS})
  INCLUDE_DIRECTORIES(${${DIRECT_DEP}_INCLUDE_DIRS})
  LINK_DIRECTORIES(${${DIRECT_DEP}_LIBRARY_DIRS})
ENDFOREACH(DIRECT_DEP ${DIRECT_DEPS})

# Set headers
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

# Set sources
SET(SOURCES
  RBGen_Params.cpp
  )

# Set addition headers and source based on package-specific 
# configuration options
IF(${PROJECT_NAME}_ENABLE_EPETRA)
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
ENDIF(${PROJECT_NAME}_ENABLE_EPETRA)

# Target library setup
ADD_LIBRARY(rbgen ${HEADERS} ${SOURCES})
LIST(APPEND LIBRARIES rbgen)
LIST(APPEND LIB_DEPS teuchos)
TARGET_LINK_LIBRARIES(rbgen ${LIB_DEPS})

# Library dependencies based on package-specific
# configuration options
IF(${PROJECT_NAME}_ENABLE_EPETRA)
  LIST(APPEND LIB_DEPS epetra)
  TARGET_LINK_LIBRARIES(rbgen ${LIB_DEPS})
ENDIF(${PROJECT_NAME}_ENABLE_EPETRA)

INSTALL(
  TARGETS rbgen
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
  COMPONENT ${PROJECT_NAME}
	)

INSTALL(
  FILES ${HEADERS}
  DESTINATION ${TRILINOS_INSTALL_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Define current package dependency variables
GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${INCLUDE_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} include directories")
GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)
SET(${PROJECT_NAME}_LIBRARY_DIRS ${${PROJECT_NAME}_LIBRARY_DIRS}  ${LIBRARY_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${LIBRARIES} CACHE INTERNAL "${PROJECT_NAME} libraries")
