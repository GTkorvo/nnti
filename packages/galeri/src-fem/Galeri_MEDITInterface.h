// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_MEDIT
#define GALERI_MEDIT

#include <iostream>
#include <iomanip>
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Galeri_AbstractGrid.h"

namespace Galeri {
namespace FiniteElements {

class MEDITInterface
{

public:

  MEDITInterface(Epetra_Comm& Comm) :
    Comm_(Comm)
  {}

  ~MEDITInterface() {}

  const Epetra_Comm& Comm() const
  {
    return(Comm_);
  }

  void Write(const AbstractGrid& data, const string& BaseName,
             const Epetra_MultiVector& Field)
  {
    ///int zzz = data.NumMyVertices();
    vector<double> coord(3);
    vector<int>    vertices(data.NumVerticesPerElement());

    const Epetra_Map& RowMap = data.RowMap();
    const Epetra_Map& VertexMap = data.VertexMap();

    Epetra_MultiVector Coord(VertexMap, 3);

    for (int i = 0 ; i < data.NumMyVertices() ; ++i)
    {
      data.VertexCoord(i, &coord[0]);
      Coord[0][i] = coord[0];
      Coord[1][i] = coord[1];
      Coord[2][i] = coord[2];
    }

    int n = 0;
    if (Field.Comm().MyPID() == 0)
      n = RowMap.NumGlobalElements();
    Epetra_Map SingleProcMap(-1, n, 0, Field.Comm());
    Epetra_MultiVector SingleProcCoord(SingleProcMap, 3);
    Epetra_MultiVector SingleProcField(SingleProcMap, 1);

    Epetra_Import CoordImporter(SingleProcMap, VertexMap);
    Epetra_Import FieldImporter(SingleProcMap, RowMap);
    SingleProcCoord.Import(Coord, CoordImporter, Insert);
    SingleProcField.Import(Field, FieldImporter, Insert);

    if (Comm().MyPID() == 0)
    {
      string FileName = BaseName + ".mesh";
      std::ofstream medit;

      medit.open(FileName.c_str());
      medit << "MeshVersionFormatted 1" << endl;
      medit << "Dimension 3" << endl;
      medit << "# mesh from ML/Trilinos" << endl << endl;
      medit << "Vertices " << data.NumGlobalVertices() << endl;

      switch (data.NumDimensions()) {
      case 2:
        for (int i = 0 ; i < SingleProcCoord.MyLength() ; ++i) {
          medit << setw(12) << setiosflags(ios::showpoint) 
            << setw(12) << SingleProcCoord[0][i] << " "
            << setw(12) << SingleProcCoord[1][i] << " 0.0 1" << endl;
        }
        break;
      case 3:
        for (int i = 0 ; i < SingleProcCoord.MyLength() ; ++i) {
          medit << setw(12) << setiosflags(ios::showpoint) 
            << setw(12) << SingleProcCoord[0][i] << " "
            << setw(12) << SingleProcCoord[1][i] << " "
            << setw(12) << SingleProcCoord[2][i] << " 1" << endl;
        }
        break;
      default:
        throw(-1);
      }

      medit.close();
    }
    Comm().Barrier();

    for (int ProcID = 0 ; ProcID < Comm().NumProc() ; ++ProcID) {

      if (Comm().MyPID() == ProcID) {

        string FileName = BaseName + ".mesh";
        std::ofstream medit(FileName.c_str(),ios::app);

        if (ProcID == 0) {
          string type = data.ElementType();

          if (type == "GALERI_TRIANGLE")
            medit << "Triangles " << data.NumGlobalElements() << endl;
          else if (type == "GALERI_QUAD")
            medit << "Quadrilaterals " << data.NumGlobalElements() << endl;
          else if (type == "GALERI_TET")
            medit << "Tetrahedra " << data.NumGlobalElements() << endl;
          else if (type == "GALERI_HEX")
            medit << "Hexahedra " << data.NumGlobalElements() << endl;
          else
          {
            cerr << "Incorrect element type (" << type << ")" << endl;
            throw(-1);
          }
        }

        for (int i = 0 ; i < data.NumMyElements() ; ++i) {
          data.ElementVertices(i, &vertices[0]);
          for (int j = 0 ; j < data.NumVerticesPerElement() ; ++j)
            medit << data.VertexMap().GID(vertices[j]) + 1 << " ";
          medit << Comm().MyPID() << endl;
        }

        if (ProcID == Comm().NumProc() - 1) 
          medit << endl << "End" << endl;

        medit.close();
      }

      Comm().Barrier();

    } // for Procs, write elements

    // ======== //
    // .bb file //
    // ======== //

    if (Comm().MyPID() == 0) {

      string BBName = BaseName + ".bb";    
      std::ofstream bb;

        bb.open(BBName.c_str());
        bb << "3 1 " << data.NumGlobalVertices() << " 2" << endl;

      for (int i = 0 ; i < SingleProcField.MyLength() ; ++i)
        bb << setiosflags(ios::showpoint) << SingleProcField[0][i] << endl;

      bb.close();

    }

  }

private:
  const Epetra_Comm& Comm_;

};

} // namespace FiniteElements
} // namespace Galeri
#endif
