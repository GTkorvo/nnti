
SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Core  core  PS  REQUIRED
  EpetraAdapters  adapters/epetra  PS  OPTIONAL
  EpetraExtAdapters  adapters/epetraext  PS  OPTIONAL
  TpetraAdapters  adapters/tpetra  PS  OPTIONAL
  )

# NOTE: The above subpackages are automatically listed as dependencies of
# Thyra!

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
