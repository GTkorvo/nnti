INCLUDE(AssertDefined)

FUNCTION(PRINT_VAR VARIBLE_NAME)
  #ASSERT_DEFINED(VARIBLE_NAME)
  MESSAGE(STATUS "${VARIBLE_NAME}='${${VARIBLE_NAME}}'")
ENDFUNCTION()
