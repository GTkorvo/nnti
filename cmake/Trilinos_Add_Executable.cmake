# $Header$

# 2008/06/09: rabartl: Below: I removed these comment blocks from other
# CMakeLists.txt files.  Which comment is correct?

# From epetra/test/CMakeLists.txt: Syntax for adding an executable:
#   TRILINOS_ADD_EXECUTABLE(<SOURCE_DIR> [TEST] [MPITEST] [EXAMPLE] [INSTALL] SOURCES <SOURCE1> <SOURCE2> ... ARGS <ARG1> <ARG2> ...)

# From epetra/test/CMakeLists.txt: Syntax for adding a test:
#   TRILINOS_ADD_EXECUTABLE(<TEST_DIR_NAME> SOURCES <SOURCE1> <SOURCE2> ... ARGS <ARG1> <ARG2> ...)

INCLUDE(Parse_Variable_Arguments)

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

# 2008/07/09: rabartl: ToDo::
#
# (*) Change the name of the current TRILINOS_ADD_EXECUTABLE(...) to
#     TRILINOS_ADD_EXECUTABLE_AND_TEST(...).
#
# (*) break these macros each into their separate *.cmake files
#


# 2008/07/09: rabartl: ToDo::
#
# (*) Change the name of the current TRILINOS_ADD_TARGET(...) to 
#     TRILINOS_ADD_EXECUTABLE(...)
#
# (*) Add an optional DIRECTORY argument and put the executable in that directory
#
# (*) Add same COMM logic for deciding whether to add an executable depending
#     on serial or COMM
#
# (*) ??? 
#
MACRO (TRILINOS_ADD_TARGET EXECUTABLE_DIR)
   PARSE_ARGUMENTS(
     PARSE  #prefix
     "SOURCES;NAME;DIRECTORY" #lists
     "INSTALL" #options
     ${ARGN} )
  
  SET (EXE_SOURCES)
  SET(EXE_BINARY_NAME "${PROJECT_NAME}-${EXEC_TYPE}-${PARSE_NAME}.exe")
  
  FOREACH( ARG ${PARSE_SOURCES} )
    SET (EXE_SOURCES ${EXE_SOURCES} ${EXECUTABLE_DIR}/${ARG})
  ENDFOREACH( ARG )                
  ADD_EXECUTABLE(${EXE_BINARY_NAME} ${EXE_SOURCES})
  IF(PARSE_DIRECTORY)
    SET_TARGET_PROPERTIES( ${EXE_BINARY_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PARSE_DIRECTORY} )
  ENDIF()
ENDMACRO()


# 2008/07/09: rabartl: ToDo:
#
# (*) Support optional DIRECTORY argument
#
# (*) Support multiple ARGS keywords
#
# (*) Support an optional POSTFIX argument for naming tests with different
# ARGS keywords

#
FUNCTION (TRILINOS_ADD_TEST EXECUTABLE_DIR)
  PARSE_ARGUMENTS(
     PARSE   #prefix
     "COMM;NAME;ARGS;PASS_REGULAR_EXPRESSION;HOST;XHOST;FAIL_REGULAR_EXPRESSION"  #lists
     "DESEND_INTO_DIR"   #options
     ${ARGN} )

  MESSAGE("PARSE_NAME = ${PARSE_NAME}")
  MESSAGE("PARSE_COMM = ${PARSE_COMM}")

  SET(DUMMY_LIST "one" "two" "three")
  MESSAGE("DUMMY_LIST = ${DUMMY_LIST}")

  # Don't add the test if the host or xhost tests don't pass

  IF(NOT PARSE_XHOST)
    SET (PARSE_XHOST NONE)
  ENDIF()    
  LIST (FIND PARSE_XHOST ${TRILINOS_HOSTNAME} INDEX_OF_HOSTNAME_IN_XHOST_LIST)           
  IF (NOT ${INDEX_OF_HOSTNAME_IN_XHOST_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  IF(NOT PARSE_HOST)
    SET (PARSE_HOST ${TRILINOS_HOSTNAME})
  ENDIF()  
  LIST (FIND PARSE_HOST ${TRILINOS_HOSTNAME} INDEX_OF_HOSTNAME_IN_HOST_LIST)                 
  IF (${INDEX_OF_HOSTNAME_IN_HOST_LIST} EQUAL -1)
    RETURN()
  ENDIF()

  # 2008/07/09: rabartl: ToDo: Above, change the logic to allow HOST and XHOST
  # to be lists!


  SET(EXE_BINARY_NAME "${PROJECT_NAME}-${EXEC_TYPE}-${PARSE_NAME}.exe")
  SET(TEST_NAME "${PROJECT_NAME}-${EXEC_TYPE}-${PARSE_NAME}")
  GET_TARGET_PROPERTY(EXECUTABLE_PATH ${EXE_BINARY_NAME} LOCATION)
  MESSAGE("mypath: ${EXECUTABLE_PATH}")
  SET(ADDED_THE_TEST OFF)    
  IF(TRILINOS_ENABLE_MPI)
    MESSAGE("MPI Mode!")
    IF (NOT PARSE_COMM)
      # If no COMM is given assume we will add the test
      SET(DO_MPI_INDEX 0)
    ELSE()
      # Else, if COMM is defined we have to find 'mpi'
      LIST (FIND PARSE_COMM "mpi" DO_MPI_INDEX)
    ENDIF()
    IF(NOT ${DO_MPI_INDEX} EQUAL -1)
      SET(TEST_NAME "${PROJECT_NAME}-${EXEC_TYPE}-${PARSE_NAME}-MPI")
      
      SET(COUNTER 0)
      FOREACH(PARSE_ARG ${PARSE_ARGS})
        #MESSAGE("ADD TEST WITH ARGS: ${PARSE_ARG}")
       
        
        
        INCREMENT(COUNTER ${COUNTER})
      ENDFOREACH()
      ADD_TEST(${TEST_NAME} ${MPI_EXECUTABLE} ${MPI_EXECUTABLE_FLAGS} ${EXECUTABLE_PATH} ${PARSE_ARGS})
      SET(ADDED_THE_TEST ON)    
    ENDIF()
  ELSE()
    MESSAGE("SERIAL Mode!")
    IF (NOT PARSE_COMM)
      # If no COMM is given assume we will add the test
      SET(DO_SERIAL_INDEX 0)
    ELSE()
      # Else, if COMM is defined we have to find 'serial'
      LIST (FIND PARSE_COMM "serial" DO_SERIAL_INDEX)
    ENDIF()
    MESSAGE("DO_SERIAL_INDEX = ${DO_SERIAL_INDEX}")
    IF(NOT ${DO_SERIAL_INDEX} EQUAL -1)
      #ADD_TEST(${TEST_NAME} ${PARSE_DIRECTORY}/${EXE_BINARY_NAME} ${PARSE_ARGS})
      ADD_TEST(${TEST_NAME} ${EXECUTABLE_PATH} ${PARSE_ARGS})
     
      SET(ADDED_THE_TEST ON)    
    ENDIF()
  ENDIF()

  # 2008/07/09: rabartl: ToDo: Above, create a macho called
  # ???ITEM_EXITS_IN_LIST??(...) to simplify logic!
    
  IF (PARSE_PASS_REGULAR_EXPRESSION AND ADDED_THE_TEST)
   SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION
     ${PARSE_PASS_REGULAR_EXPRESSION})
  ENDIF()
  
  IF (PARSE_FAIL_REGULAR_EXPRESSION AND ADDED_THE_TEST)
   SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES FAIL_REGULAR_EXPRESSION
     ${PARSE_FAIL_REGULAR_EXPRESSION})
  ENDIF()

ENDFUNCTION(TRILINOS_ADD_TEST)


MACRO(INCREMENT var_name input)
  # Increment the input variable (must be in [0,8]) and store the result in var_name.
  SET(${var_name} ${increment${input}})
  IF(NOT DEFINED ${var_name})
    MESSAGE(FATAL_ERROR "Could not increment. Input ${input} out of range 0-8?")
  ENDIF(NOT DEFINED ${var_name})
ENDMACRO(INCREMENT)


SET(increment0 1)
SET(increment1 2)
SET(increment2 3)
SET(increment3 4)
SET(increment4 5)
SET(increment5 6)
SET(increment6 7)
SET(increment7 8)
SET(increment8 9)
