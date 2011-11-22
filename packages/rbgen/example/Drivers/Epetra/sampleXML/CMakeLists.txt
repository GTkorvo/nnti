
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopySampleXMLExampleFiles
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/rbgen/example/Drivers/Epetra/sampleXML
  SOURCE_FILES sampleInput.xml README 
  )
