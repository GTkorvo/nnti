
INCLUDE(ParseVariableArguments)
INCLUDE(PrintVar)


SET(CMAKE_EXECUTABLE_SUFFIX ".exe")


MACRO(PACKAGE_PRIVATE_ADD_TEST_SET_PASS_PROPERTY TEST_NAME_IN)

  IF (PARSE_PASS_REGULAR_EXPRESSION)
    SET_TESTS_PROPERTIES(${TEST_NAME_IN} PROPERTIES PASS_REGULAR_EXPRESSION
      ${PARSE_PASS_REGULAR_EXPRESSION})
  ENDIF()

  IF (PARSE_FAIL_REGULAR_EXPRESSION)
    SET_TESTS_PROPERTIES(${TEST_NAME_IN} PROPERTIES FAIL_REGULAR_EXPRESSION
      ${PARSE_FAIL_REGULAR_EXPRESSION})
  ENDIF()

  IF (PARSE_STANDARD_PASS_OUTPUT)
    SET_TESTS_PROPERTIES(${TEST_NAME_IN} PROPERTIES PASS_REGULAR_EXPRESSION
      "End Result: TEST PASSED")
  ENDIF()

ENDMACRO()


#
# Function that converts a complete string of command-line arguments
# into a form that ADD_TEST(...) can correctly deal with.
#
# The main thing this function does is to replace spaces ' ' with
# array separators ';' since this is how ADD_TEST(...) expects to deal
# with command-line arguments, as array arguments.  However, this
# function will not do a replacement of ' ' with ';' if a quote is
# active.  This allows you to pass in quoted arguments and have them
# treated as a single argument.
#

FUNCTION(CONVERT_CMND_ARG_STRING_TO_ADD_TEST_ARG_ARRAY CMND_ARG_STRING ARG_ARRAY_VARNAME)
  
  #MESSAGE("CONVERT_CMND_ARG_STRING_TO_ADD_TEST_ARG_ARRAY")
  #PRINT_VAR(CMND_ARG_STRING)
  #PRINT_VAR(ARG_ARRAY_VARNAME)

  STRING(LENGTH ${CMND_ARG_STRING} STR_LEN)
  #PRINT_VAR(STR_LEN)

  MATH(EXPR STR_LAST_IDX "${STR_LEN}-1")

  SET(NEWSTR)

  SET(ACTIVE_QUOTE OFF)

  FOREACH(IDX RANGE ${STR_LAST_IDX})

    STRING(SUBSTRING ${CMND_ARG_STRING} ${IDX} 1 STR_CHAR)
    #PRINT_VAR(STR_CHAR)

    IF (STR_CHAR STREQUAL "\"")
      IF (NOT ACTIVE_QUOTE)
        SET(ACTIVE_QUOTE ON)
      ELSE()
        SET(ACTIVE_QUOTE OFF)
      ENDIF()
      #PRINT_VAR(ACTIVE_QUOTE)
    ENDIF()

    IF (NOT STR_CHAR STREQUAL " ")
      SET(NEWSTR "${NEWSTR}${STR_CHAR}")
    ELSE()
      IF (ACTIVE_QUOTE)
        SET(NEWSTR "${NEWSTR}${STR_CHAR}")
      ELSE()
        SET(NEWSTR "${NEWSTR};")
      ENDIF()
    ENDIF()

  ENDFOREACH()
   
  #PRINT_VAR(NEWSTR)

  SET(${ARG_ARRAY_VARNAME} ${NEWSTR} PARENT_SCOPE)

ENDFUNCTION()