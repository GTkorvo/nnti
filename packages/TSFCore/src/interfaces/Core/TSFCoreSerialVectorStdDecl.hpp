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

// /////////////////////////////////////////////////////////////////
// TSFCoreSerialVectorDecl.hpp

#ifndef TSFCORE_SERIAL_VECTOR_STD_DECL_HPP
#define TSFCORE_SERIAL_VECTOR_STD_DECL_HPP

#include "TSFCoreSerialVectorBaseDecl.hpp"

namespace TSFCore {

///
/** General extendable implementation of serial vectors.
 *
 * This class can be used either as a view of a vector data or as a
 * storage for vector data (with any underlying storage type).
 *
 * To create with storage with the dimension of <tt>dim</tt> just call
 * the constructor <tt>SerialVectorStd(dim)</tt> or after construction
 * you can call <tt>this->initialize(dim)</tt>.
 *
 * To simply create a view of a vector <tt>v</tt> with stride
 * <tt>vs</tt>, without ownership just call
 * <tt>SerialVectorStd(Teuchos::rcp(v,false),vs)</tt> or after
 * construction call
 * <tt>this->initialize(Teuchos::rcp(v,false),vs)</tt>.
 *
 * To use another storage type, such as an
 * <tt>std::vector<Scalar></tt>, construct as:
 *
 \code

 template<class Scalar>
 Teuchos::RefCountPtr<Vector<Scalar> > STLVectorSpace<Scalar>::createMember() const
 {
   Teuchos::RefCountPtr<std::vector<Scalar> > stl_v = Teuchos::rcp( new std::vector<Scalar>(dim_) );
   Teuchos::RefCountPtr<Scalar> v = Teuchos::rcp(&(*stl_v)[0],false);
	 Teuchos::set_extra_data( stl_v, "stl::vector", &v );
	 return Teuchos::rcp( new SerialVectorStd<Scalar>( v, 1, dim_, Teuchos::rcp(this,false) ) );
 }

 \endcode
 *
 * \ingroup TSFCore_adapters_serial_concrete_std_grp
 */
template<class Scalar>
class SerialVectorStd : public SerialVectorBase<Scalar> {
public:

	/** @name Constructors/initializers */
	//@{

	/** \brief Calls <tt>this->initialize(vecSpc)</tt>.
	 */
	SerialVectorStd(
		const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
		);
	/** \brief Calls <tt>this->initialize(dim)</tt>.
	 */
	SerialVectorStd(
		const Index dim = 0
		);
	/** \brief Calls <tt>this->initialize(v,vs,dim,vecSpc)</tt>.
	 */
	SerialVectorStd(
		const Teuchos::RefCountPtr<Scalar>                      &v
		,const Index                                            vs
		,const Index                                            dim
		,const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc = Teuchos::null
		);
	/** \brief Call <tt>this->initialize(v,vs,vecSpc)</tt> with internally dynamically allocated data <tt>v</tt>.
	 */
	void initialize(
		const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc
		);
	/** \brief Call <tt>this->initialize(v,vs,true)</tt> with internally dynamically allocated data <tt>v</tt>.
	 */
	void initialize(
		const Index dim
		);
	/** \brief Initialize with storage.
	 *
	 * @param  v      [in] Smart pointer to array of storage that <tt>*this</tt> will represent.
	 * @param  vs     [in] Stride for the storage in <tt>v[]</tt> (see Postconditions).
	 * @param  dim    [in] Number of elements in <tt>v[]</tt> this this will represent (see Postconditions).
	 * @param  vecSpc
	 *                [in] Smart pointer to a <tt>VectorSpace</tt> object that will be used to represent the
	 *                vector space for <tt>*this</tt>.  If <tt>vecSpc.get()==NULL</tt> on input, then
	 *                a <tt>SerialVectorSpace</tt> object of dimension <tt>dim</tt> is allocated for this
	 *                role.  The default is <tt>Teuchos::null</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>dim == vecSpc->dim()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>vecSpc->createMember()</tt> must create vectors that are compatible
	 *      with <tt>*this</tt> (i.e. <tt>getSubVector()</tt>, <tt>commitSubVector()</tt> behave the same as with
	 *      this class).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> [<tt>vecSpc.get()!=NULL</tt>] <tt>vecSpc.get() == this->space().get()</tt>
	 * <li> [<tt>vecSpc.get()==NULL</tt>] <tt>dynamic_cast<const SerialVectorSpaceStd<Scalar>*>(this->space().get()) != NULL</tt>
	 * <li> <tt>this->space()->dim() == dim</tt>
	 * <li> <tt>this->getRCPtr().get() == v.get()</tt>
	 * <li> <tt>this->getPtr() == v.get()</tt>
	 * <li> <tt>this->getStride() == vs</tt>
	 * </ul>
	 *
	 * Note that this function is declared virtual so that subclasses
	 * can override it to be informed whenever <tt>*this</tt> vector
	 * is resized.  An override should call this function as
	 * <tt>this->SerialVectorStd<Scalar>::initialize(...)</tt>.
	 */
	virtual void initialize(
		const Teuchos::RefCountPtr<Scalar>                      &v
		,const Index                                            vs
		,const Index                                            dim
		,const Teuchos::RefCountPtr<const VectorSpace<Scalar> > &vecSpc = Teuchos::null
		);

	//@}

	/** @name Accessors (inlined for minimal overhead) */
	//@{

	///
	Teuchos::RefCountPtr<Scalar> getRCPtr();
	///
	Teuchos::RefCountPtr<const Scalar> getRCPtr() const;
	///
	Scalar* getPtr();
	///
	const Scalar* getPtr() const;
	///
	Index getStride() const;
	///
	Index getDim() const;
	
	//@}

	/** @name Overridden from SerialVectorBase */
	//@{
	///
	void getData( Scalar** values, Index* stride );
	///
	void commitData( Scalar** values );
	///
	void getData( const Scalar** values, Index* stride ) const;
	///
	void freeData( const Scalar** values ) const;
	//@}

	/** @name Overridden from Vector */
	//@{
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > space() const;
	//@}

  /** @name Overridden from Teuchos::Describable */
  //@{
	///
	std::string description() const;
	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	Teuchos::RefCountPtr<Scalar>                            v_;
	Index                                                   vs_;
	Index                                                   dim_;
	Teuchos::RefCountPtr<const VectorSpace<Scalar> >        space_serial_;

	// ////////////////////////////////
	// Private member functions

	void free_mem();

	// Not defined and not to be called
	SerialVectorStd(const SerialVectorStd&);
	SerialVectorStd& operator=(const SerialVectorStd&);

}; // end class SerialVectorStd

// /////////////////////////////////////////////////////
// Inline members

template<class Scalar>
inline
Teuchos::RefCountPtr<Scalar> SerialVectorStd<Scalar>::getRCPtr()
{
	return v_;
}

template<class Scalar>
inline
Teuchos::RefCountPtr<const Scalar> SerialVectorStd<Scalar>::getRCPtr() const
{
	return v_;
}

template<class Scalar>
inline
Scalar* SerialVectorStd<Scalar>::getPtr()
{
	return v_.get();
}

template<class Scalar>
inline
const Scalar* SerialVectorStd<Scalar>::getPtr() const
{
	return v_.get();
}

template<class Scalar>
inline
Index SerialVectorStd<Scalar>::getStride() const
{
	return vs_;
}	

template<class Scalar>
inline
Index SerialVectorStd<Scalar>::getDim() const
{
	return dim_;
}	

} // end namespace TSFCore

#endif // TSFCORE_SERIAL_VECTOR_STD_DECL_HPP
