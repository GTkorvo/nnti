//@HEADER
/*
************************************************************************

              EpetraExt: Extended Linear Algebra Services Package 
                Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef EPETRAEXT_BLOCKMULTIVECTOR_H
#define EPETRAEXT_BLOCKMULTIVECTOR_H

#include "Epetra_MultiVector.h" 

#include <vector>

//! EpetraExt::BlockMultiVector: A class for constructing a distributed block multivector

/*! The EpetraExt::BlockMultiVector allows construction of a block multivector made up of Epetra_MultiVector blocks as well as access to the full systems as a Epetra_MultiVector.  It derives from and extends the Epetra_MultiVector class

<b>Constructing EpetraExt::BlockMultiVector objects</b>

*/    

namespace EpetraExt {

class BlockMultiVector: public Epetra_MultiVector {
 public:

  //@{ \name Constructors/Destructor.
  //! BlockMultiVector constuctor with one block row per processor.
  /*! Creates a BlockMultiVector object and allocates storage.  
    
	\param In
	BaseMap - Map determining local structure, can be distrib. over subset of proc.'s
	\param In 
	GlobalMap - Full map describing the overall global structure, generally generated by the construction of a BlockCrsMatrix object
	\param In 
	NumVectors - Number of vectors in object
  */
  BlockMultiVector( const Epetra_BlockMap & BaseMap, const Epetra_BlockMap & GlobalMap,
 int NumVectors );
  
  //! BlockMultiVector constuctor
  /*! Creates a BlockMultiVector object and allocates storage.  
    
	\param In
	BaseMap - Map determining local structure, can be distrib. over subset of proc.'s
	\param In 
	GlobalMap - Full map describing the overall global structure, generally generated by the construction of a BlockCrsMatrix object
	\param In 
	NumBlocks - Number of blocks on this processor
	\param In 
	NumVectors - Number of vectors in object
  */
  BlockMultiVector( const Epetra_BlockMap & BaseMap, const Epetra_BlockMap & GlobalMap,
                    int NumBlocks, int NumVectors );

  //! Copy constructor.
  BlockMultiVector( const BlockMultiVector & MV );

  //! Destructor
  virtual ~BlockMultiVector();
  //@}
  
  //! Block Access
  Epetra_MultiVector & Block( int i = 0 ) { return *(Blocks_[i]); }
  const Epetra_MultiVector & Block( int i = 0 ) const { return *(Blocks_[i]); }
	
 protected:

  void AllocateBlocks_();
  void DeleteBlocks_();

  Epetra_BlockMap BaseMap_;

  std::vector<Epetra_MultiVector*> Blocks_;

  int NumBlocks_;

  std::vector< double** > Ptrs_;

};

} //namespace EpetraExt

#endif /* EPETRAEXT_BLOCKMULTIVECTOR_H */
