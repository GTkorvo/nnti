/* @HEADER@ */
/* @HEADER@ */

#ifndef TSF_PARAMUTILS_H
#define TSF_PARAMUTILS_H

#include "TSFConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"


namespace TSFExtended
{
Teuchos::ParameterList mergeParamLists(
  const Teuchos::ParameterList& pDefault, 
  const Teuchos::ParameterList& pIn);
}

#endif
