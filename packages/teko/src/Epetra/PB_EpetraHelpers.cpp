#include "PB_EpetraHelpers.hpp"

// Thyra Includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

// Epetra includes
#include "Epetra_Vector.h"

// EpetraExt includes
#include "EpetraExt_ProductOperator.h"
#include "EpetraExt_MatrixMatrix.h"

// PB includes
#include "PB_EpetraOperatorWrapper.hpp"
#include "PB_Utilities.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;
using Teuchos::null;

namespace PB {
namespace Epetra {

/** \brief Convert an Epetra_Vector into a diagonal linear operator.
  *
  * Convert an Epetra_Vector into a diagonal linear operator. 
  *
  * \param[in] ev  Epetra_Vector to use as the diagonal
  * \param[in] map Map related to the Epetra_Vector
  * \param[in] lbl String to easily label the operator
  *
  * \returns A diagonal linear operator using the vector
  */
const Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraDiagOp(const RCP<const Epetra_Vector> & ev,const Epetra_Map & map,
                                                                   const std::string & lbl)
{
   const RCP<const Thyra::VectorBase<double> > thyraVec  // need a Thyra::VectorBase object
         = Thyra::create_Vector(ev,Thyra::create_VectorSpace(rcpFromRef(map)));
   Teuchos::RCP<Thyra::LinearOpBase<double> > op 
         = Teuchos::rcp(new Thyra::DefaultDiagonalLinearOp<double>(thyraVec));
   op->setObjectLabel(lbl);
   return op;
}

/** \brief Fill a Thyra vector with the contents of an epetra vector. This prevents the
  *
  * Fill a Thyra vector with the contents of an epetra vector. This prevents the need
  * to reallocate memory using a create_MultiVector routine. It also allows an aritrary
  * Thyra vector to be filled.
  *
  * \param[in,out] spmdMV Multi-vector to be filled.
  * \param[in]     mv     Epetra multi-vector to be used in filling the Thyra vector.
  */    
void fillDefaultSpmdMultiVector(Teuchos::RCP<Thyra::DefaultSpmdMultiVector<double> > & spmdMV,
                                Teuchos::RCP<Epetra_MultiVector> & epetraMV)
{
   // first get desired range and domain
   const RCP<const Thyra::SpmdVectorSpaceBase<double> > range  = spmdMV->spmdSpace();
   const RCP<const Thyra::ScalarProdVectorSpaceBase<double> > domain 
         = rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<double> >(spmdMV->domain());

   TEUCHOS_ASSERT(domain->dim()==epetraMV->NumVectors());

   // New local view of raw data
   double *localValues; int leadingDim;
   if(epetraMV->ConstantStride() )
      epetraMV->ExtractView( &localValues, &leadingDim );
   else
      TEST_FOR_EXCEPT(true); // ToDo: Implement views of non-contiguous mult-vectors!

   // Build the MultiVector
   spmdMV->initialize(range, domain,
                     Teuchos::arcp(localValues,0,leadingDim*epetraMV->NumVectors(),false),
                     leadingDim);

   // make sure the Epetra_MultiVector doesn't disappear prematurely
   Teuchos::set_extra_data<RCP<Epetra_MultiVector> >(epetraMV,"Epetra_MultiVector",Teuchos::outArg(spmdMV));
}

/** \brief Build a vector of the dirchlet row indicies. 
  *
  * Build a vector of the dirchlet row indicies. That is, record the global
  * index of any row that is all zeros except for $1$ on the diagonal.
  *
  * \param[in]     rowMap   Map specifying which global indicies this process examines 
  * \param[in]     mat      Matrix to be examined
  * \param[in,out] indicies Output list of indicies corresponding to dirchlet rows.
  */
void identityRowIndicies(const Epetra_Map & rowMap, const Epetra_CrsMatrix & mat,std::vector<int> & outIndicies)
{
   int maxSz = mat.GlobalMaxNumEntries();
   double values[maxSz];
   int indicies[maxSz];

   // loop over elements owned by this processor
   for(int i=0;i<rowMap.NumMyElements();i++) {
      bool rowIsIdentity = true;
      int sz = 0;
      int rowGID = rowMap.GID(i);
      mat.ExtractGlobalRowCopy(rowGID,maxSz,sz,values,indicies);

      // loop over the columns of this row
      for(int j=0;j<sz;j++) {
         int colGID = indicies[j];

         // look at row entries
         if(colGID==rowGID) rowIsIdentity &= values[j]==1.0;
         else               rowIsIdentity &= values[j]==0.0;

         // not a dirchlet row...quit
         if(not rowIsIdentity)
            break;
      }

      // save a row that is dirchlet
      if(rowIsIdentity)
         outIndicies.push_back(rowGID);
   }
}

/** \brief Zero out the value of a vector on the specified
  *        set of global indicies.
  *
  * Zero out the value of a vector on the specified set of global
  * indicies. The indicies here are assumed to belong to the calling
  * process (i.e. zeroIndicies $\in$ mv.Map()).
  *
  * \param[in,out] mv           Vector whose entries will be zeroed
  * \param[in]     zeroIndicies Indicies local to this process that need to be zeroed
  */
void zeroMultiVectorRowIndicies(Epetra_MultiVector & mv,const std::vector<int> & zeroIndicies)
{
   int colCnt = mv.NumVectors();
   std::vector<int>::const_iterator itr;
 
   // loop over the indicies to zero
   for(itr=zeroIndicies.begin();itr!=zeroIndicies.end();++itr) {
 
      // loop over columns
      for(int j=0;j<colCnt;j++)
         mv.ReplaceGlobalValue(*itr,j,0.0);
   }
}

/** \brief Constructor for a ZeroedOperator.
  *
  * Build a ZeroedOperator based on a particular Epetra_Operator and
  * a set of indicies to zero out. These indicies must be local to this
  * processor as specified by RowMap().
  *
  * \param[in] zeroIndicies Set of indices to zero out (must be local).
  * \param[in] op           Underlying epetra operator to use.
  */
ZeroedOperator::ZeroedOperator(const std::vector<int> & zeroIndicies,
                               const Teuchos::RCP<const Epetra_Operator> & op)
   : zeroIndicies_(zeroIndicies), epetraOp_(op)
{ }

//! Perform a matrix-vector product with certain rows zeroed out
int ZeroedOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
   int result = epetraOp_->Apply(X,Y);

   // zero a few of the rows
   zeroMultiVectorRowIndicies(Y,zeroIndicies_);

   return result;
}

} // end namespace Epetra
} // end namespace PB
