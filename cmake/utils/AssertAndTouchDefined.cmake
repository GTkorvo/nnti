FUNCTION(ASSERT_AND_TOUCH_DEFINED VARS)
  FOREACH(VAR ${VARS})
    IF(NOT DEFINED ${VAR})
      MESSAGE(SEND_ERROR "Error, the variable ${VAR} is not defined!")
    ELSE()
      # Read the varaible so that it will register as being read!
      SET(DUMMY_VAR ${${VAR}})
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()
