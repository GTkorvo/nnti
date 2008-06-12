# $Header$

ADD_SUBDIRECTORY(Lena)
ADD_SUBDIRECTORY(Moon)
ADD_SUBDIRECTORY(SampleXML)

# Setup include directories and library dependencies
INCLUDE_DIRECTORIES(${${PROJECT_NAME}_INCLUDE_DIRS})
LINK_LIBRARIES(${${PROJECT_NAME}_LIBRARIES})

#ADD_EXECUTABLE(RBGenDriver_EpetraMV ${${PROJECT_NAME}_SOURCE_DIR}/src/RBGenDriver_EpetraMV.cpp)

#TARGET_LINK_LIBRARIES(RBGenDriver_EpetraMV rbgen teuchos epetra)

