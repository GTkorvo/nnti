//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include <EDT_CrsGraph_View.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_BlockMap.h>

#include <vector>

using std::vector;

namespace EpetraExt {

CrsGraph_View::
~CrsGraph_View()
{
  if( newObj_ ) delete newObj_;
}

CrsGraph_View::NewTypeRef
CrsGraph_View::
operator()( OriginalTypeRef orig )
{
  origObj_ = &orig;

  //Error, must be local indices
  assert( !orig.IndicesAreGlobal() );

  //test maps, new map must be left subset of old
  const Epetra_BlockMap & oRowMap = orig.RowMap();
  const Epetra_BlockMap & oColMap = orig.ColMap();

  int oNumRows = oRowMap.NumMyElements();
  int oNumCols = oRowMap.NumMyElements();
  int nNumRows = NewRowMap_->NumMyElements();
  int nNumCols = 0;
  if( NewColMap_ ) nNumCols = NewColMap_->NumMyElements();

  bool matched = true;
  for( int i = 0; i < nNumRows; ++i )
    matched = matched && ( oRowMap.GID(i) == NewRowMap_->GID(i) );
  if( nNumCols )
    for( int i = 0; i < nNumCols; ++i )
      matched = matched && ( oColMap.GID(i) == NewColMap_->GID(i) );

  if( !matched ) cout << "EDT_CrsGraph_View: Bad Row or Col Mapping\n";
  assert( matched );

  //intial construction of graph
  vector<int> numIndices( nNumRows );
  vector<int*> indices( nNumRows );
  for( int i = 0; i < nNumRows; ++i )
  {
    orig.ExtractMyRowView( i, numIndices[i], indices[i] );
    int j = 0;
    if( nNumCols )
    {
      while( j < numIndices[i] && NewColMap_->GID(indices[i][j]) != -1 ) ++j;
      numIndices[i] = j;
    }
  }

  Epetra_CrsGraph * newGraph( new Epetra_CrsGraph( View,
                                                   *NewRowMap_,
                                                   *NewColMap_,
                                                   &numIndices[0] ) );

  //insert views of row indices
  for( int i = 0; i < nNumRows; ++i )
    newGraph->InsertMyIndices( i, numIndices[i], indices[i] );

  newGraph->TransformToLocal();

  newObj_ = newGraph;

  return *newGraph;
}

} // namespace EpetraExt

