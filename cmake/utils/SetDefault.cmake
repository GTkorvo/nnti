MACRO(SET_DEFAULT VAR)
  IF ("${${VAR}}" STREQUAL "")
    SET(${VAR} ${ARGN})
  ENDIF()
ENDMACRO()
