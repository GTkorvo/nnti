#include "Epetra/PB_StridedMappingStrategy.hpp"
#include "Epetra/PB_InterlacedEpetra.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

namespace PB {
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
StridedMappingStrategy::StridedMappingStrategy(const std::vector<int> & vars,const RCP<const Epetra_Map> & map,
                                               const Epetra_Comm & comm)
{
   rangeMap_ = map;
   domainMap_ = map;
   buildBlockTransferData(vars, rangeMap_,comm);
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
void StridedMappingStrategy::copyEpetraIntoThyra(const Epetra_MultiVector& X,
                                                 const Teuchos::Ptr<Thyra::MultiVectorBase<double> > & thyra_X,
                                                 const PB::Epetra::EpetraOperatorWrapper & eow) const
{
   int count = X.NumVectors(); 

   std::vector<RCP<Epetra_MultiVector> > subX;
   std::vector<RCP<const Epetra_MultiVector> > subY;

   // allocate vectors to copy into
   PB::Epetra::buildSubVectors(blockMaps_,subX,count);

   // copy source vector to X vector
   PB::Epetra::one2many(subX,X,blockImport_);

   // convert subX to an array of multi vectors
   Teuchos::Array<RCP<const Thyra::MultiVectorBase<double> > > thyra_subX;
   for(int i=0;i<blockMaps_.size();i++)
      thyra_subX.push_back(Thyra::create_MultiVector(subX[i],Thyra::productVectorSpaceBase(eow.getThyraOp()->domain())->getBlock(i)));
  
   // build product multivector
   const RCP<const Thyra::DefaultProductVectorSpace<double> > pvs 
         = rcp_dynamic_cast<const Thyra::DefaultProductVectorSpace<double> >(eow.getThyraOp()->domain());
   Teuchos::ptr_dynamic_cast<Thyra::DefaultProductMultiVector<double> >(thyra_X)->initialize(pvs,thyra_subX);
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
void StridedMappingStrategy::copyThyraIntoEpetra(const RCP<const Thyra::MultiVectorBase<double> > & thyra_Y,
                                                 Epetra_MultiVector& Y,
                                                 const PB::Epetra::EpetraOperatorWrapper & eow) const
{
   std::vector<RCP<const Epetra_MultiVector> > subY;
   RCP<const Thyra::DefaultProductMultiVector<double> > prod_Y 
         = rcp_dynamic_cast<const Thyra::DefaultProductMultiVector<double> >(thyra_Y);

   int numBlocks = prod_Y->productSpace()->numBlocks();

   // convert thyra product vector to subY
   for(int i=0;i<blockMaps_.size();i++)
      subY.push_back(Thyra::get_Epetra_MultiVector(*blockMaps_[i].second,prod_Y->getMultiVectorBlock(i)));

   // copy solution vectors to Y vector
   PB::Epetra::many2one(Y,subY,blockExport_);
}

// this is the core routine that builds the maps
// and importers/exporters neccessary for all the
// transfers. Currently it simply calls out to the
// interlaced epetra functions. (Comment: this
// routine should probably be private or protected
// ... it is basically the meat of the constructor)
//
//    arguments:
//       vars - Vector describing the blocking of variables
//       baseMap - basic map to use in the transfers
//       comm    - Epetra_Comm object
//
void StridedMappingStrategy::buildBlockTransferData(const std::vector<int> & vars,const Teuchos::RCP<const Epetra_Map> & baseMap, const Epetra_Comm & comm)
{
   int numGlobals = baseMap->NumGlobalElements();

   // build maps and exporters/importers
   PB::Epetra::buildSubMaps(numGlobals,vars,comm,blockMaps_);
   PB::Epetra::buildExportImport(*baseMap, blockMaps_, blockExport_,blockImport_);
}

// Builds a blocked Thyra operator that uses the strided
// mapping strategy to define sub blocks.
//
//    arguments:
//       mat - Epetra_CrsMatrix with FillComplete called, this
//             matrix is assumed to be square, with the same
//             range and domain maps
//    returns: Blocked Thyra linear operator with sub blocks
//             defined by this mapping strategy
//
const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
StridedMappingStrategy::buildBlockedThyraOp(const RCP<const Epetra_CrsMatrix> & crsContent,const std::string & label) const
{
   int dim = blockMaps_.size();

   RCP<Thyra::DefaultBlockedLinearOp<double> > A = Thyra::defaultBlockedLinearOp<double>();

   A->beginBlockFill(dim,dim);
   for(int i=0;i<dim;i++) {
      for(int j=0;j<dim;j++) {
         // label block correctly
         std::stringstream ss;
         ss << label << "_" << i << "," << j;

         // build the blocks and place it the right location
         A->setBlock(i,j,Thyra::epetraLinearOp(PB::Epetra::buildSubBlock(i,j,*crsContent,blockMaps_),ss.str()));
      }
   } // end for i
   A->endBlockFill();

   return A;
}

} // end namespace Epetra
} // end namespace PB
