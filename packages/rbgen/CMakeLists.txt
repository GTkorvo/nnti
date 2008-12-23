
INCLUDE(PackageMacros)

PACKAGE(RBGen)

SET(${PACKAGE_NAME}_VERSION "0.1d")

ADD_SUBDIRECTORY(src)

# There are not tests and examples yet!

PACKAGE_POSTPROCESS()
