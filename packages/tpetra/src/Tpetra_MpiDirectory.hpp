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

#ifndef TPETRA_MPIDIRECTORY_HPP
#define TPETRA_MPIDIRECTORY_HPP

#include "Tpetra_Directory.hpp"
#include "Tpetra_Object.hpp"
#include "Tpetra_ElementSpace.hpp"

namespace Tpetra {

  //! Tpetra::MpiDirectory: This class is an MPI implementation of Directory. Its interface allows ElementSpace and BlockElementSpace objects to reference non-local elements.
  
  /*! For ElementSpace objects, a Directory object must be created by a call to
    the Comm createDirectory method.  The Directory is needed to allow referencing
    of non-local elements.
		For BlockElementSpace objects, a Directory should be created and used through the 
		ElementSpace accessor.
    
		This class currently has two constructors, one that takes an ElementSpace object, and a copy constructor.
  */
  
  template<typename OrdinalType>
  class MpiDirectory : public Object, public virtual Directory<OrdinalType> {
  public:
    
    //@{ \name Constructors/Destructor.
    //! constructor
    MpiDirectory(ElementSpace<OrdinalType> const& elementSpace)	
      : Object("Tpetra::Directory[MPI]") 
      , ElementSpace_(&elementSpace) 
    {};
    
    //! copy constructor
    MpiDirectory(MpiDirectory<OrdinalType> const& directory)
      : Object(directory.label()) 
      , ElementSpace_(directory.ElementSpace_) 
    {};
    
    //! destructor.
    ~MpiDirectory() {};
    //@}
    
    //@{ \name Query method.
    //! getDirectoryEntries : Returns image and local id info for non-local ElementSpace entries
    /*! Given a list of Global Entry IDs, this function returns the list of
      image IDs and local IDs on the owning memory image that correspond
      to the list of entries.  If LocalEntries is 0, then local IDs are 
      not returned.  If EntrySizes is nonzero, it will contain a list of corresponding 
      element sizes for the requested global entries.
      \param In
      numEntries - Number of Global IDs being passed in.
      \param In
      globalEntries - List of Global IDs being passed in.
      \param InOut
      images - User allocated array of length at least NumEntries.  On return contains list of images
      owning the Global IDs in question.
      \param InOut
      localEntries - User allocated array of length at least NumEntries.  On return contains the local ID of
      the global on the owning image. If LocalEntries is zero, no local ID information is returned.
    */
    void getDirectoryEntries(OrdinalType numEntries, OrdinalType const* globalEntries, 
                             OrdinalType* images, OrdinalType* localEntries) const 
    {
      throw reportError("This method is not implemented yet.", -1);
    };
    //@}
    
  private:
    ElementSpace<OrdinalType> const* ElementSpace_;
    
    //! Assignment operator (declared but not defined, do not use)
    MpiDirectory<OrdinalType>& operator = (MpiDirectory<OrdinalType> const& Source);
    
  }; // class MpiDirectory
  
} // namespace Tpetra

#endif // TPETRA_MPIDIRECTORY_HPP
