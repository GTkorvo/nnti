// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef _TPETRA_SERIALDISTRIBUTOR_HPP_
#define _TPETRA_SERIALDISTRIBUTOR_HPP_

#include "Tpetra_Object.hpp"
#include "Tpetra_Distributor.hpp"

namespace Tpetra {

//! Tpetra::SerialDistributor:  The Tpetra Serial implementation of the Tpetra::Distributor Gather/Scatter Setup Class.
/*! The SerialDistributor class is an Serial implement of Tpetra::Distributor that is essentially a trivial class
    since a serial machine is a trivial parallel machine.
		An SerialDistributor object is actually produced by calling a method in the Tpetra::SerialPlatform class.

		Most SerialDistributor methods throw an error of -1, since they should never be called.
*/

	template<typename PacketType, typename OrdinalType>
	class SerialDistributor : public Object, public virtual Distributor<PacketType, OrdinalType> {
  public:

  //@{ \name Constructor/Destructor

  //! Default Constructor.
  SerialDistributor() : Object("Tpetra::Distributor[Serial]") {};

  //! Copy Constructor
  SerialDistributor(SerialDistributor<PacketType, OrdinalType> const& plan) : Object(plan.label()) {};

  //! Clone method
	Distributor<PacketType, OrdinalType>* clone() {
		Distributor<PacketType, OrdinalType>* distributor = static_cast<Distributor<PacketType, OrdinalType>*>
			(new SerialDistributor<PacketType, OrdinalType>(*this)); 
		return(distributor); 
	};

  //! Destructor.
  virtual ~SerialDistributor() {};
  //@}

	//@{ \name Gather/Scatter Constructors
  //! Create Distributor object using list of Image IDs to send to
  void createFromSends(OrdinalType const& numExportIDs, OrdinalType const* exportImageIDs,
											 bool const& deterministic, OrdinalType& numRemoteIDs ) 
		{throw reportError("This method should never be called.", -1);};
	//! Create Distributor object using list of Image IDs to receive from
	void createFromRecvs(OrdinalType const& numRemoteIDs, OrdinalType const* remoteGIDs, 
											 OrdinalType const* remoteImageIDs, bool const& deterministic, 
											 OrdinalType& numExportIDs, OrdinalType*& exportGIDs, 
											 OrdinalType*& exportImageIDs)
		{throw reportError("This method should never be called.", -1);};
	//@}
	
	//@{ \name Constant Size
	//! do
  void doPostsAndWaits       (PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs) 
		{throw reportError("This method should never be called.", -1);};
	//! doReverse
  void doReversePostsAndWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};

	//! doPosts
  void doPosts(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//! doWaits
  void doWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};

	//! doReversePosts
  void doReversePosts(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//! doReverseWaits
  void doReverseWaits(PacketType* export_objs, OrdinalType const& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//@}

	//@{ \name Non-Constant Size
	//! do
  void doPostsAndWaits       (PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//! doReverse
  void doReversePostsAndWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};

	//! doPosts
  void doPosts(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//! doWaits
  void doWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};

	//! doReversePosts
  void doReversePosts(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//! doReverseWaits
  void doReverseWaits(PacketType* export_objs, OrdinalType const*& obj_size, PacketType* import_objs)
		{throw reportError("This method should never be called.", -1);};
	//@}

	//@{ \name I/O Methods
	//! print method inherited from Object
  void print(ostream& os) const {os << label();};
	//! printInfo method inherited from Distributor
  void printInfo(ostream& os) const {print(os);};
	//@}

}; // class SerialDistributor

} // namespace Tpetra

#endif // _TPETRA_SERIALDISTRIBUTOR_HPP_
