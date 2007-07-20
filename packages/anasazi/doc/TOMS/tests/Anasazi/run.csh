#!/bin/csh 

echo "*** Running ./BKS_2DQ1.exe --nx=$1 --maxRestarts=3 --tol=0"
echo ""
time ./BKS_2DQ1.exe --nx=$1 --maxRestarts=3 --tol=0 --echo-command-line
