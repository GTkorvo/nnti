PROJECT(RBGen)

SET(RBGEN_VERSION "0.1d")

OPTION(RBGEN_ENABLE_EPETRA
  "Enable RBGen features that use epetra."
  ON
  )

#OPTION(RBGEN_ENABLE_TESTS
#  "Build RBGen tests."
#  OFF
#  )

#OPTION(RBGEN_ENABLE_EXAMPLES
#  "Build RBGen examples."
#  OFF
#  )

ADD_SUBDIRECTORY(src)

# Note ... currently no tests are defined, so this is a no-op
IF(RBGEN_ENABLE_TESTS)
  ADD_SUBDIRECTORY(test)
ENDIF(RBGEN_ENABLE_TESTS)

# Note ... currently no examples are defined, so this is a no-op
IF(RBGEN_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(example)
ENDIF(RBGEN_ENABLE_EXAMPLES)

