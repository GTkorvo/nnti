
INCLUDE(TribitsPackageMacros)

TRIBITS_PACKAGE(RBGen)

SET(${PACKAGE_NAME}_VERSION "0.1d")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(example)

# There are not tests and examples yet!

TRIBITS_PACKAGE_POSTPROCESS()
