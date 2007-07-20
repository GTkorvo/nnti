#!/bin/csh 

set args = " --nx=$1 --ncv=$2 --maxRestarts=0 --tol=0 --nonherm"
echo "*** Running ./BKS_2DQ1.exe $args"
echo ""
time ./BKS_2DQ1.exe $args --nonherm --echo-command-line
