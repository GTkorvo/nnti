// //////////////////////////////////////////////////////////////////////
// SparsePtrElement.h

#ifndef SPARSE_PTR_ELEMENT_H
#define SPARSE_PTR_ELEMENT_H

#include "SparseLinAlgPackTypes.h"

namespace SparseLinAlgPack {

///
/** Sparse pointer element type.
  *
  * This class abstracts a sparse element of a templated
  * type.  It is ment to be used in a sparse vector.  It
  * has a pointer to the value of the element.
  *
  * The default assignment operator and copy constructor
  * are allowed.
  */
template <class T_Indice, class T_Value>
class SparsePtrElement {
public:
	/** @name Public Typedefs. */
	//@{

	///
	typedef T_Value							value_type;
	///
	typedef T_Indice						indice_type;

	//@}

	/** @name Constructors */
	//@{

	/// Construct uninitialized (poiner to value set to zero) (#indice() == 0#).
	SparsePtrElement() : indice_(0), pvalue_(0)
	{}

	/// Construct with a pointer to the value and indice set
	SparsePtrElement(indice_type indice, value_type* pvalue) : indice_(indice), pvalue_(pvalue)
	{}
	
	//@}

	/** @name Value and indice access */
	//@{ 

	///
	value_type& value()
	{
		return *pvalue_;
	}
	///
	value_type value() const
	{
		return *pvalue_;
	}
	///
	indice_type indice() const
	{
		return indice_;
	}
	/// Change the indice
	void change_indice(indice_type indice)
	{
		indice_ = indice;
	}
	/// Change the element pointer
	void change_value_ptr(value_type* pvalue)
	{
		pvalue_ = pvalue;
	}

	//@}
private:
	indice_type				indice_;
	value_type*				pvalue_;

};	// end class SparsePtrElement

} // end namespace SparseLinAlgPack 

#endif // SPARSE_PTR_ELEMENT_H
