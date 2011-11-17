
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyLenaExampleFiles
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/rbgen/example/Drivers/Epetra/Lena
  SOURCE_FILES lenaInputAnasazi.xml lenaInputIncSVD.xml lenaInputStSVD.xml lenaInput.xml lena.mtx README
  )
