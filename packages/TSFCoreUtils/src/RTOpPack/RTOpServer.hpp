// ////////////////////////////////////////////////////
// RTOpServer.hpp

#ifndef RTOP_SERVER_HPP
#define RTOP_SERVER_HPP

#include "RTOpServerDecl.hpp"

namespace RTOpPack {

template<class Scalar>
void RTOpServer<Scalar>::add_op_factory(
	const MemMngPack::ref_count_ptr<MemMngPack::AbstractFactory<RTOpPack::RTOpT<Scalar> > > &op_factory
	)
{
	// ToDo: RAB: 20030620: Validate op_factory properly!
	op_factories_[op_factory->create()->op_name()] = op_factory;
}

template<class Scalar>
MemMngPack::ref_count_ptr<MemMngPack::AbstractFactory<RTOpPack::RTOpT<Scalar> > >
RTOpServer<Scalar>::get_op_factory( const char op_name[] ) const
{
	typename op_factories_t::const_iterator itr = op_factories_.find(op_name);
	THROW_EXCEPTION(
		itr == op_factories_.end(), std::logic_error
		,"RTOpServer<Scalar>::get_op_factory(...): Error, an operator factory with the "
		"operator name \'" << op_name << "\' does not exist!"
		);
	return itr->second;
}

template<class Scalar>
void RTOpServer<Scalar>::print_op_factories(std::ostream& o) const
{
	using std::setw;
	const int w = 30;
	o << "\nRTOpServer<Scalar>::print_op_factories(...): RTOp operator factories currently registered\n\n" << std::left;
	o << setw(w) << "Operator name" << "Operator type" << std::endl;
	o << setw(w) << "-------------" << "-------------" << std::endl;
	for( typename op_factories_t::const_iterator itr = op_factories_.begin(); itr != op_factories_.end(); ++itr ) {
		o << setw(w) << itr->first << typeid(*itr->second->create()).name() << std::endl;
	}
	o << std::endl;
}

} // namespace RTOpPack

#endif // RTOP_SERVER_HPP
