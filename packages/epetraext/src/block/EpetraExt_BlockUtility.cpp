//@HEADER
// ************************************************************************
// 
//               EpetraExt: Extended Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "EpetraExt_BlockUtility.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

namespace EpetraExt {

using std::vector;

//==============================================================================
Epetra_CrsGraph * BlockUtility::GenerateBlockGraph(
        const Epetra_CrsGraph & BaseGraph,
        const vector< vector<int> > & RowStencil,
        const vector<int> & RowIndices,
        const Epetra_Comm & GlobalComm ) 
{
  const Epetra_BlockMap & BaseMap = BaseGraph.RowMap();
  int MaxGID = BaseMap.MaxAllGID();
  int BaseIndex = BaseMap.IndexBase();

  //Setup block offset
  int Offset = 1;
  while( Offset < MaxGID ) Offset *= 10;

  //Get Base Global IDs
  int NumBlockRows = RowIndices.size();
  int Size = BaseMap.NumMyElements();
  int TotalSize = NumBlockRows * Size;
  vector<int> GIDs(Size);
  BaseMap.MyGlobalElements( &GIDs[0] );

  vector<int> GlobalGIDs( TotalSize );
  for( int i = 0; i < NumBlockRows; ++i )
  {
    for( int j = 0; j < Size; ++j )
      GIDs[i*Size+j] += RowIndices[i] * Offset;
  }

  int GlobalSize;
  GlobalComm.SumAll( &TotalSize, &GlobalSize, 1 );

  Epetra_Map GlobalMap( GlobalSize, TotalSize, &GlobalGIDs[0], BaseIndex, GlobalComm );

  int MaxIndices = BaseGraph.MaxNumIndices();
  vector<int> Indices(MaxIndices);
  int NumIndices;

  Epetra_CrsGraph * GlobalGraph = new Epetra_CrsGraph( Copy, 
                               dynamic_cast<Epetra_BlockMap&>(GlobalMap),
                               0 );

  for( int i = 0; i < NumBlockRows; ++i )
  {
    int StencilSize = RowStencil[i].size();
    for( int j = 0; j < Size; ++j )
    {
      int BaseRow = BaseMap.GID(j);
      int GlobalRow = GlobalMap.GID(j+i*Size);

      BaseGraph.ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
      for( int k = 0; k < StencilSize; ++k )
      {
        int ColOffset = (RowIndices[i]+RowStencil[i][k]) * Offset;
        if( k > 0 ) ColOffset -= (RowIndices[i]+RowStencil[i][k-1]) * Offset;

        for( int l = 0; l < NumIndices; ++l )
          Indices[l] += ColOffset;

        GlobalGraph->InsertGlobalIndices( GlobalRow, NumIndices, &Indices[0] );
      }
    }
  }

  GlobalGraph->TransformToLocal();
  
  return GlobalGraph;
}

} //namespace EpetraExt
