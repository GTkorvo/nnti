

FUNCTION(JOIN  OUTPUT_STRING_VAR  SEP_STR)
  SET(OUTPUT_STRING "")
  FOREACH(STRING_VAL ${ARGN})
    IF (OUTPUT_STRING STREQUAL "")
      SET(OUTPUT_STRING "${STRING_VAL}")
    ELSE()
      SET(OUTPUT_STRING "${OUTPUT_STRING}${SEP_STR}${STRING_VAL}")
    ENDIF()
  ENDFOREACH()
  SET(${OUTPUT_STRING_VAR} "${OUTPUT_STRING}" PARENT_SCOPE)
ENDFUNCTION()
