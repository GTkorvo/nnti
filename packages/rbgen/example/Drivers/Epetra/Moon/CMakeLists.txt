
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyMoonExampleFiles
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/rbgen/example/Drivers/Epetra/Moon
  SOURCE_FILES moonInput.xml moon.mtx README
  )
