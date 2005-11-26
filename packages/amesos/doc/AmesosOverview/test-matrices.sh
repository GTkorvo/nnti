#!/bin/bash

EXE=$HOME/Trilinos/G4_MPI/packages/amesos/example/example_AmesosFactory_HB.exe
DIR=$HOME/projects/MATRICES/FIDAP
OUTPUT=res-superlu
SOLVER=Amesos_Superlu

/bin/rm -f $OUTPUT

for i in $DIR/*.rua
do
  $EXE $i $SOLVER 2>&1 >> $OUTPUT
done
