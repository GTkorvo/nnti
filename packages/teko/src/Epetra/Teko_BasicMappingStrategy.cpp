#include "Epetra/Teko_BasicMappingStrategy.hpp"
#include "Epetra/Teko_EpetraHelpers.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

namespace Teko {
namespace Epetra {

// Creates a strided mapping strategy. This class is useful
// for breaking up nodally ordered matrices (i.e. the unknowns
// in a FEM problem are ordered [u0,v0,p0,u1,v1,p1,...]). Current
// implimentation only supports a fixed number of variables
//
//    arguments: 
//       vars - Number of different variables 
//       map  - original Epetra_Map to be broken up
//       comm - Epetra_Comm object related to the map
//
BasicMappingStrategy::BasicMappingStrategy(const Teuchos::RCP<const Epetra_Map> & rMap,
                                           const Teuchos::RCP<const Epetra_Map> & dMap, const Epetra_Comm & comm)
{
   rangeMap_ = rMap;
   domainMap_ = dMap;
}

// Virtual function defined in MappingStrategy.  This copies
// an Epetra_MultiVector into a Thyra::MultiVectorBase with
// blocking handled by the strides defined in the constructor.
//
//   arguments:
//      X       - source Epetra_MultiVector
//      thyra_X - destination Thyra::MultiVectorBase
//      eow     - Operator that defines the transition...this may
//                be removed in the future
//
void BasicMappingStrategy::copyEpetraIntoThyra(const Epetra_MultiVector& X,
                                               const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & thyra_X,
                                               const Teko::Epetra::EpetraOperatorWrapper & eow) const
{
   // perform a simple copy
   RCP<Thyra::DefaultSpmdMultiVector<double> > vec 
            = rcp_dynamic_cast<Thyra::DefaultSpmdMultiVector<double> >(Teuchos::rcpFromRef(*thyra_X)); 
   Teuchos::RCP<Epetra_MultiVector> ptrX = Teuchos::rcp_const_cast<Epetra_MultiVector>(Teuchos::rcpFromRef(X));
   fillDefaultSpmdMultiVector(vec,ptrX);
}

// Virtual function defined in MappingStrategy.  This copies
// an Epetra_MultiVector into a Thyra::MultiVectorBase with
// blocking handled by the strides defined in the constructor.
//
//   arguments:
//      thyra_Y - source Thyra::MultiVectorBase
//      Y       - destination Epetra_MultiVector
//      eow     - Operator that defines the transition...this may
//                be removed in the future
//
void BasicMappingStrategy::copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> > & thyra_Y,
                                                 Epetra_MultiVector& Y,
                                                 const Teko::Epetra::EpetraOperatorWrapper & eow) const
{
   RCP<const Epetra_MultiVector> eSrc = Thyra::get_Epetra_MultiVector(*rangeMap(),thyra_Y);

   Y = *eSrc;
}

} // end namespace Epetra
} // end namespace Teko
