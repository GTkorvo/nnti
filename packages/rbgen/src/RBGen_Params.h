#ifndef RBGEN_CREATE_PARAMS_H
#define RBGEN_CREATE_PARAMS_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace RBGen {

Teuchos::RefCountPtr<Teuchos::ParameterList> createParams( int argc, char* argv[] );

}
#endif
