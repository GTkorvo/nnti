
INCLUDE(TrilinosPackageHelpers)

TRILINOS_PACKAGE(RBGen)

SET(${PACKAGE_NAME}_VERSION "0.1d")

MESSAGE("RBGen_ENABLE_Epetra = ${RBGen_ENABLE_Epetra}") 

ADD_SUBDIRECTORY(src)

# There are not tests and examples yet!

TRILINOS_PACKAGE_POSTPROCESS()
