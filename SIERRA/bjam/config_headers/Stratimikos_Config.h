
#include "TrilinosPlatform_config.h"

#ifndef STRATIMIKOS_FUNC_TIME_MONITOR
#  define STRATIMIKOS_TEUCHOS_TIME_MONITOR
#  define STRATIMIKOS_FUNC_TIME_MONITOR(FUNCNAME) \
     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, STRATIMIKOS)
#  define STRATIMIKOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF) \
     TEUCHOS_FUNC_TIME_MONITOR_DIFF(FUNCNAME, DIFF)
#endif

