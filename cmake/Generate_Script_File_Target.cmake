
#
# Generate a script file given a script body string and a set of optional
# file dependencies
#

MACRO(GENERATE_SCRIPT_FILE_TARGET SCRIPT_FILE_NAME SCRIPT_BODY)
  
  SET(TEST_SCRIPT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${SCRIPT_FILE_NAME}")

  FILE(WRITE ${TEST_SCRIPT_FILE}.tmp "${SCRIPT_BODY}")

  ADD_CUSTOM_COMMAND(
    OUTPUT ${SCRIPT_FILE_NAME}
    DEPENDS ${ARGN}
    COMMAND cp ARGS ${SCRIPT_FILE_NAME}.tmp ${SCRIPT_FILE_NAME}
    COMMAND chmod ARGS a+x ${SCRIPT_FILE_NAME}
    COMMAND rm ARGS ${SCRIPT_FILE_NAME}.tmp
    )

  ADD_CUSTOM_TARGET( ${SCRIPT_FILE_NAME}-target ALL
    DEPENDS ${SCRIPT_FILE_NAME} )

ENDMACRO()
