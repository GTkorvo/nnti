// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"

#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Triangle.h"
#include "Galeri_grid_Quad.h"
#include "Galeri_grid_Loadable.h"
#include "Galeri_viz_MEDIT.h"

using namespace Galeri;

// ============================================================================
// When run on two processors, the code creates the following grid:
//
// [2] -- [5] -- [8]
//  |    / |    / |
//  |   /  |   /  |
//  |  /   |  /   |
//  | /    | /    |
// [1] -- [4] -- [7]
//  |      |      |
//  |      |      |
//  |      |      |
//  |      |      |
// [0] -- [3] -- [6]
//
// which is composed by 2 quadrilateral elements, and 4 triangular elements.
// The quads have size h x h. h is the minimal length of triangle sides as
// well. When run with more processors, additional elements are added
// on the right of the above picture, so that the total number of quads is
// comm.NumProc(), and that of triangles 2 * comm.NumProc().
//
// Then, one each processor we define a vector, living the of grid vertices,
// and we produce a VTK output file. This file can be visualized, for example,
// with MaYaVi.
//
// \author Marzio Sala, ETH
//
// \date Last modified on Aug-06.
// ============================================================================

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Galeri::core::Workspace::setNumDimensions(2);

  // domain1 contains the quads, domain2 the triangles
  grid::Loadable domain1(comm, -1, 1, "Quad");
  grid::Loadable domain2(comm, -1, 2, "Triangle");

  int vertexOffset = 3 * comm.MyPID();
  int elementOffset = 3 * comm.MyPID();
  double h = 1.0;
  double xOffset = h * comm.MyPID();

  // domain1, start with connectivity.
  domain1.setGlobalConnectivity(elementOffset, 0, vertexOffset + 0);
  domain1.setGlobalConnectivity(elementOffset, 1, vertexOffset + 3);
  domain1.setGlobalConnectivity(elementOffset, 2, vertexOffset + 4);
  domain1.setGlobalConnectivity(elementOffset, 3, vertexOffset + 1);

  domain1.freezeConnectivity();

  // x-coordinates
  domain1.setGlobalCoordinates(vertexOffset + 0, 0, xOffset);
  domain1.setGlobalCoordinates(vertexOffset + 3, 0, xOffset + h);
  domain1.setGlobalCoordinates(vertexOffset + 4, 0, xOffset + h);
  domain1.setGlobalCoordinates(vertexOffset + 1, 0, xOffset);

  // y-coordinates
  domain1.setGlobalCoordinates(vertexOffset + 0, 1, 0.0);
  domain1.setGlobalCoordinates(vertexOffset + 3, 1, 0.0);
  domain1.setGlobalCoordinates(vertexOffset + 4, 1, h);
  domain1.setGlobalCoordinates(vertexOffset + 1, 1, h);

  // now domain2, start with connectivity
  domain2.setGlobalConnectivity(elementOffset, 0, vertexOffset + 1);
  domain2.setGlobalConnectivity(elementOffset, 1, vertexOffset + 4);
  domain2.setGlobalConnectivity(elementOffset, 2, vertexOffset + 5);

  domain2.setGlobalConnectivity(elementOffset + 1, 0, vertexOffset + 1);
  domain2.setGlobalConnectivity(elementOffset + 1, 1, vertexOffset + 5);
  domain2.setGlobalConnectivity(elementOffset + 1, 2, vertexOffset + 2);

  domain2.freezeConnectivity();

  // x-coordinates
  domain2.setGlobalCoordinates(vertexOffset + 1, 0, xOffset);
  domain2.setGlobalCoordinates(vertexOffset + 4, 0, xOffset + h);
  domain2.setGlobalCoordinates(vertexOffset + 5, 0, xOffset + h);
  domain2.setGlobalCoordinates(vertexOffset + 2, 0, xOffset);

  // y-coordinates
  domain2.setGlobalCoordinates(vertexOffset + 1, 1, h);
  domain2.setGlobalCoordinates(vertexOffset + 4, 1, h);
  domain2.setGlobalCoordinates(vertexOffset + 5, 1, 2 * h);
  domain2.setGlobalCoordinates(vertexOffset + 2, 1, 2 * h);

  Epetra_Vector vector(domain1.getVertexMap());
  vector.Random();

  Galeri::viz::MEDIT::write(domain1, "domain1", vector);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
