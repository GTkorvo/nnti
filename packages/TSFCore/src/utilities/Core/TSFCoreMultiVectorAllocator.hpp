// ///////////////////////////////////////////////////////////////
// TSFCoreMultiVectorAllocator.hpp

#ifndef TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP
#define TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP

#include "TSFCoreVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {

///
/** Allocator class to be used with <tt>MemMngPack::AbstractFactoryStd</tt> to create
 * <tt>MultiVector</tt> objects of a given size.
 */
template<class Scalar>
class MultiVectorAllocator {
public:
	///
	MultiVectorAllocator() : numMembers_(0) {}
	///
	typedef Teuchos::RefCountPtr<MultiVector<Scalar> >  ptr_t;         // required!
	///
	MultiVectorAllocator( const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vs, int numMembers )
		: vs_(vs), numMembers_(numMembers)
		{
#ifdef _DEBUG
			TEST_FOR_EXCEPTION( vs.get()==NULL, std::logic_error, "Error!" );
#endif			
		}
	///
	const ptr_t allocate() const { return vs_->createMembers(numMembers_); }  // required!
private:
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >      vs_;
	int                                                        numMembers_;
};

} // namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP
