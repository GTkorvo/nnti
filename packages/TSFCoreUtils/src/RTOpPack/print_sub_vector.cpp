// //////////////////////////////////////////////////////////////////
// print_sub_vector.cpp

#include "print_sub_vector.hpp"

std::ostream& RTOpPack::output(
	std::ostream& o, const SubVector& v
	,bool print_dim , bool newline
	)
{
	int w = o.width(0) - 1; // get the set width (minus 1 since a space is inserted)
	if( print_dim )
//		o << std::setw(0) << std::left << v.subDim() << std::endl << std::right;
 		o << std::setiosflags(std::ios::left) << std::setw(0) << v.subDim() 
		  << std::endl << std::setiosflags(std::ios::right);
	// RAB: 20030916: ToDo: Fix the above by hacking std::left and std::right in config header!
	const RTOp_value_type  *v_val        = v.values();
	const ptrdiff_t        v_val_s       = v.stride();
	for( RTOp_index_type i = 1; i <= v.subDim(); ++i, v_val+=v_val_s ) {
		// insert a space to be sure there is white space
		// inbetween adjacent elements.
		o << " " << std::setw(w) << (*v_val) << ":" << i + v.globalOffset();
	}
	if(newline) o << std::endl;
	return o;
}
