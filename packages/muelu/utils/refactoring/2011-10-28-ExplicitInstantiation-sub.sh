python /home/jngaida/dev/MueLu/src/preCopyrightTrilinos/muelu/utils/refactoring/2011-10-28-ExplicitInstantiation.py -i $1 -o `dirname $1`/`basename $1 .hpp`_decl.hpp -l $2
#tkdiff $1 `dirname $1`/`basename $1 .hpp`_decl.hpp