#!/bin/bash

# Returns the list of files to be copied over

for file_base in `grep 'input' apdx_TeuchosMemMngSummary.tex | perl -pe 's/\\\\input{(.+)}/$1.*/'` ; do echo "$file_base.*" ; done
for file_base in `grep 'includegraphics' apdx_TeuchosMemMngSummary.tex | perl -pe 's/\\\\includegraphics.+{(.+)}.*/$1/'` ; do echo "$file_base.*" ; done
