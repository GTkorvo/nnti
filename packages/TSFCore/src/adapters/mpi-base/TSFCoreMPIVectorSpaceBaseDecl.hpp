// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIVectorSpaceBaseDecl.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_BASE_DECL_HPP
#define TSFCORE_MPI_VECTOR_SPACE_BASE_DECL_HPP

#include "TSFCoreVectorSpaceStdBase.hpp"

namespace TSFCore {

///
/** <tt>%VectorSpace</tt> node subclass for all MPI-based vectors with
 * contiguous local storage.
 *
 * This interface collaborates with the classes <tt>MPIVetorBase</tt>
 * and <tt>MPIMultiVetorBase</tt> to implement MPI-based parallel (or
 * serial of course) vectors and multi-vectors that allows different
 * concrete implementations to mix together.
 *
 * Specifically, these classes are designed to handle three different
 * use case:
 * <ul>
 * <li> Serial vectors on one processor with <tt>localSubDim()==dim()</tt>.
 *      In this case, the <tt>MPI_Comm</tt> returned from <tt>mpiComm()</tt> can
 *      be <tt>MPI_COMM_NULL</tt>.
 * <li> Distributed parallel vectors on two or more processors with
 *      <tt>localSubDim() < dim()</tt>.  This implementation assumes
 *      that all of the elements are stored contiguously on each
 *      processor and that there are no ghost elements as described
 *      below.
 * <li> Locally replicated vectors on one or more processor.  This
 *      case is similar to serial vectors except that we have
 *      <tt>localSubDim()==dim()</tt> even if there is more than
 *      one processor.  This case is checked for so that the reduction
 *      operations are performed correctly.
 * </ul>
 *
 * This interface provides all the information necessary to implement
 * <tt>MPIVectorBase::applyOp()</tt> in all of the above described use
 * cases.  This interface returns an MPI communicator (of which all
 * compatible vector spaces must have the same communicator obviously)
 * through the method <tt>mpiComm()</tt>.
 *
 * <A NAME="MPIVectorSpaceBase_Vector_layout"></A>
 * <b>%Vector data layout:</b>
 *
 * For the case of a distributed parallel vector, this interface base
 * class assumes that vector data is partitioned to processors in
 * contiguous chunks of dense subvectors.  To spell this out, let
 * <tt>v</tt> be the local vector that is sorted on this processor and
 * let <tt>g</tt> be the global vector.  Then these two vectors are
 * related (using one-based indexing) as:

 \verbatim
    v(k) == g(k + this->localOffset()), for k = 1...this->localSubDim()

 \endverbatim

 * Any type of mapping of vector data to processors that can not be
 * interpreted in this way can not rely on this base class for
 * interoperability.  Note that as long as the elements in a processor
 * are partitioned to unique processors and no ghost elements are
 * present, the actual indexes used by the application with these
 * vectors is immaterial.  The indexes associated with this set of
 * interfaces, however, are only meaningful to abstract numerial
 * algorithms and provide an arbitrary label for certain types of
 * coordinate-dependent operations (like required in an active-set
 * method for optimization).  Therefore, as long as the underlying
 * vector represents a unique partitioning of elements, these classes
 * can be used.  There is a default implementation of
 * <tt>localOffset()</tt> that automatically assumes this contiguous
 * mapping of elements to processors and in general this should not be
 * changed.
 *
 * <b>Notes to subclass developers:</b>
 *
 * The pure virtual methods <tt>mpiComm()</tt> and
 * <tt>localSubDim()</tt> defined in this interface along with the
 * pure virtual methods <tt>dim()</tt> and <tt>createMember()</tt> are
 * the only methods that must be overridden.
 *
 * If <tt>this</tt> this is in an uninitialized state then
 * <tt>localSubDim()</tt> should return <tt>0</tt>.
 *
 * If it is possible that the mapping of vector elements to processors
 * is not as described above, then the subclass should override the
 * <tt>mapCode()</tt> and <tt>isCompatible()</tt> methods as described
 * above and below but this should never be necessary or desirable to
 * do.
 *
 * If optimized implementations of multi-vectors can be supported,
 * then the <tt>createMembers()</tt> method should also be overridden.
 *
 * This class defines a very general default implementation for
 * <tt>smallVecSpcFcty()</tt> that returns a
 * <tt>MPIVectorSpaceFactoryStd</tt> object creates
 * <tt>MPIVectorSpaceStd</tt> <tt>VectorSpace</tt> objects that create
 * <tt>MPIMultiVectorStd</tt> <tt>MultiVector</tt> objects (and
 * <tt>VectorMultiVector</tt>-wrapped <tt>MPIMultiVectorStd</tt>
 * <tt>Vector</tt> objects).  This implementation should be very
 * appropriate for many different concrete implementations.
 *
 * <b>Note:</b> It is very important that subclasses call the
 * <tt>updateState()</tt> function whenever the state of
 * <tt>*this</tt> changes in a way that might affect the behavior of
 * any of the public member functions.  For example, if a different
 * value of <tt>localSubDim()</tt> will be returned the next time it
 * is called by a client, then <tt>%updateState()</tt> needs to be
 * called by the subclass.  Clients should never need to worry about
 * this function and that is why <tt>%updateState()</tt> is declared
 * protected.
 * 
 */
template<class Scalar>
class MPIVectorSpaceBase : public VectorSpaceStdBase<Scalar> {
public:

	///
	MPIVectorSpaceBase();
	
	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	///
	/** Returns the MPI communicator.
	 */
	virtual MPI_Comm mpiComm() const = 0;
	///
	/** Returns the number of local elements stored on this processor.
	 *
	 * If <tt>this</tt> this is uninitialized then <tt>localSubDim()</tt>
	 * returns <tt>0</tt>.
	 */
 	virtual Index localSubDim() const = 0;

	//@}

	/** @name Virtual methods with default implementations */
	//@{

	///
	/** Returns the offset for the local sub-vector stored on this
	 * processor.
	 *
	 * This method has a default implementation which just assigns
	 * this offset based on counting up <tt>localSubDim() on each
	 * processor and then setting <tt>localOffset()</tt> by the rank
	 * of the processor.  For example, if there are 5 elements in
	 * process 0 and 4 elements in process rank, then
	 * <tt>localOffset</tt> on each of these processors will be set
	 * as: <tt>localOffset=0</tt> on process 0, <tt>localOffset=5</tt>
	 * on process 1, <tt>localOffset=9</tt> on process 2 and so on.
	 */
	virtual Index localOffset() const;
	///
	/** Returns the code for the mapping of elements to processors.
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>this->localSubDim() > 0</tt>] <tt>this->mapCode() > 0</tt>.
	 * <li> [<tt>this->localSubDim() <= 0</tt>] <tt>this->mapCode() <= 0</tt>.
	 * </ul>
	 *
	 * This method takes the data <tt>mpiComm()</tt>, <tt>numProc</tt>
	 * (where <tt>numProc</tt> is returned from
	 * <tt>MPI_Comm_size(this->mpiComm(),&numProc)</tt>,
	 * <tt>localOffset()</tt> or <tt>localSubDim()</tt> on each
	 * processor and then uses it to compute a value for
	 * <tt>mapCode</tt> (using a single global reduction if
	 * <tt>numProc > 1</tt>) which is returned from this function.
	 *
	 * The value returned from this default implementation of this
	 * method must not be changed or this approach breaks down.  The
	 * only reason for overridding this method is for the subclass to
	 * be alerted of <em>when</em> this method is called but not
	 * <em>what</em> is returned from this method.  If a subclass
	 * developer does not understand what this means then <b>don't</b>
	 * override this method!
	 *
	 * The default implementation will always return <tt>return >
	 * 0</tt> (unless <tt>this</tt> is uninitialized) so that if this
	 * method is overriden to return <tt>return <= </tt> then this is
	 * a flag that the underlying vector map does not satisfy the
	 * assumptions of this vector space interface and vectors that are
	 * in <tt>*this</tt> vector space can not collaborate with other
	 * MPI-based vector implementations.
	 */
	virtual Index mapCode() const;

	//@}

	/** @name Overridden from VectorSpace */
	//@{

	///
	/** Returns true if all of the elements are stored on one processor.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return == (dim()==localSubDim())</tt>.
	 * </ul>
	 */
 	bool isInCore() const;

  ///
  /** Returns a <tt>MPIVectorSpaceFactoryStd</tt> object that has been given <tt>mpiComm()</tt>.
   */
	Teuchos::RefCountPtr< const VectorSpaceFactory<Scalar> > smallVecSpcFcty() const;

	///
	/** Checks the general compatibility of parallel (or serial on one
	 * processor) MPI-based vector spaces.
	 *
	 * @return Returns true if <tt>*this</tt> and <tt>vecSpace</tt>
	 * are both serial in-core vectors or if <tt>vecSpc</tt> is of
	 * type <tt>MPIVectorSpaceBase<Scalar></tt> and both <tt>*this</tt>
	 * and <tt>vecSpc</tt> have the same MPI communicators and the same
	 * mapping of elements to processors.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return = ( this->isInCore() &&
	 *      vecSpc->isInCore() ) || ( (mpiVecSpc = dynamic_cast<const
	 *      MPIVectorSpaceBase<Scalar>*>(&vecSpc)) && this->mpiComm() ==
	 *      mpiVecSpc->mpiComm() && this->mapCode() ==
	 *      mpiVecSpc->mpiComm())</tt>.
	 *
	 * </ul>
	 *
	 * If the mapping of vector elements to processors is not as
	 * described
	 * <A HREF="classTSFCore_1_1MPIVectorSpaceBase.html#MPIVectorSpaceBase_Vector_layout>above</A>
	 * then this method should be overridden in a way that is specific
	 * to the vector implementation.
	 */
 	bool isCompatible(const VectorSpace<Scalar>& vecSpc) const;
	
	//@}

protected:

	///
	/** This function must be called whenever the state of
	 * <tt>this</tt> changes and some internal state must be updated.
	 *
	 * Note that calling this function will involve one or more global
	 * reductions being called if this is parallel vector space so it
	 * should only be called when needed by subclasses.
	 *
	 * Usually, this operation only needs to be called once for every
	 * *new* parallel vector space constructed and very few parallel
	 * vector spaces will be created per application usually.
	 */
	virtual void updateState();

private:

	// //////////////////////////////////////
	// Private data members

	Index     mapCode_;    // < 0 is a flag that everything needs initialized
	bool      isInCore_;
	Index     defaultLocalOffset_;

  Teuchos::RefCountPtr< const VectorSpaceFactory<Scalar> >  smallVecSpcFcty_;
	
}; // end class MPIVectorSpaceBase

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_SPACE_BASE_DECL_HPP
