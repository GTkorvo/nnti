# $Header$

MACRO(TRILINOS_ADD_EXECUTABLE EXECUTABLE_DIR)

  SET(KEYWORD )
  SET(EXECUTABLE_NAME ${EXECUTABLE_DIR})
  SET(EXECUTABLE_SOURCES )
  SET(EXECUTABLE_ARGS )
  SET(INSTALL_EXECUTABLE 0)
  SET(TEST_EXECUTABLE 0)
  SET(MPI_TEST_EXECUTABLE 0)

  FOREACH(ARGUMENT ${ARGN})
    IF(ARGUMENT STREQUAL "NAME")
      SET(KEYWORD ${ARGUMENT})
    ELSEIF(ARGUMENT STREQUAL "SOURCES")
      SET(KEYWORD ${ARGUMENT})
    ELSEIF(ARGUMENT STREQUAL "ARGS")
      SET(KEYWORD ${ARGUMENT})
    ELSEIF(ARGUMENT STREQUAL "TYPE")
    ELSEIF(ARGUMENT STREQUAL "INSTALL")
      SET(INSTALL_EXECUTABLE 1)
    ELSEIF(ARGUMENT STREQUAL "TEST")
      SET(TEST_EXECUTABLE 1)
      SET(EXEC_TYPE "Test")
    ELSEIF(ARGUMENT STREQUAL "EXAMPLE")
      SET(EXEC_TYPE "Example")
    ELSEIF(ARGUMENT STREQUAL "MPITEST")
      SET(MPI_TEST_EXECUTABLE 1)
      SET(EXEC_TYPE "Test")
    ELSE(ARGUMENT STREQUAL "NAME")
      IF(KEYWORD STREQUAL "NAME")
        SET(EXECUTABLE_NAME ${ARGUMENT})
      ELSEIF(KEYWORD STREQUAL "SOURCES")
        SET(EXECUTABLE_SOURCES ${EXECUTABLE_SOURCES} ${EXECUTABLE_DIR}/${ARGUMENT})
      ELSEIF(KEYWORD STREQUAL "ARGS")
        SET(EXECUTABLE_ARGS ${EXECUTABLE_ARGS} ${ARGUMENT})
      ELSE(KEYWORD STREQUAL "NAME")
        MESSAGE(SEND_ERROR "TRILINOS_ADD_EXECUTABLE() - unknown keyword: ${KEYWORD}")
      ENDIF(KEYWORD STREQUAL "NAME")
    ENDIF(ARGUMENT STREQUAL "NAME")
  ENDFOREACH(ARGUMENT)

  SET(EXECUTABLE_BINARY "${PROJECT_NAME}-${EXEC_TYPE}-${EXECUTABLE_NAME}")
  SET(TEST_NAME "${PROJECT_NAME}-${EXEC_TYPE}-${EXECUTABLE_NAME}")
  SET(MPI_TEST_NAME "${PROJECT_NAME}-${EXEC_TYPE}-${EXECUTABLE_NAME}-MPI")

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_DIR})
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/${EXECUTABLE_DIR})

  ADD_EXECUTABLE(${EXECUTABLE_BINARY} ${EXECUTABLE_SOURCES})
  GET_TARGET_PROPERTY(EXECUTABLE_PATH ${EXECUTABLE_BINARY} LOCATION)

  IF(INSTALL_EXECUTABLE)
    INSTALL(TARGETS ${EXECUTABLE_BINARY} RUNTIME DESTINATION bin)
  ENDIF(INSTALL_EXECUTABLE)

  IF(TEST_EXECUTABLE)
    ADD_TEST(${TEST_NAME} ${EXECUTABLE_BINARY} ${EXECUTABLE_ARGS})
  ENDIF(TEST_EXECUTABLE)

  IF(MPI_TEST_EXECUTABLE)
    ADD_TEST(${MPI_TEST_NAME} ${MPI_EXECUTABLE} ${MPI_EXECUTABLE_FLAGS} ${EXECUTABLE_PATH} ${EXECUTABLE_ARGS})
  ENDIF(MPI_TEST_EXECUTABLE)

ENDMACRO(TRILINOS_ADD_EXECUTABLE) 

