#!/bin/bash

TEUCHOS_MM_DIR=$1
echo $TEUCHOS_MM_DIR

filesList=`././get_teuchos_mm_report_files_list.sh`
#echo $filesList

for file in $filesList ; do
  echo "Copy $file"
  cp $TEUCHOS_MM_DIR/$file* .
  done