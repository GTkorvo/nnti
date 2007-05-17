#ifndef RBGEN_CREATE_PARAMS_H
#define RBGEN_CREATE_PARAMS_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace RBGen {

  //! Create a Teuchos::ParameterList from an XML file.
  /*!
   */
  Teuchos::RefCountPtr<Teuchos::ParameterList> createParams( const std::string& filename );

  //! Extract the filename list from a Teuchos::ParameterList.
  /*!
   */
  Teuchos::RefCountPtr<std::vector<std::string> > genFileList( const Teuchos::ParameterList& params );

}
#endif
