/* @HEADER@ */
// ************************************************************************
// 
//                              Sundance
//                 Copyright (2005) Sandia Corporation
// 
// Copyright (year first published) Sandia Corporation.  Under the terms 
// of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government 
// retains certain rights in this software.
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
// Questions? Contact Kevin Long (krlong@sandia.gov), 
// Sandia National Laboratories, Livermore, California, USA
// 
// ************************************************************************
/* @HEADER@ */

#include "SundanceGeomUtils.hpp"
#include "SundanceMesh.hpp"
#include "SundancePoint.hpp"
#include "SundanceSet.hpp"
#include <queue>

using namespace Sundance;
using namespace Teuchos;
using namespace Sundance;


namespace Sundance
{
  bool cellContainsPoint(const Mesh& mesh, int cellDim, int cellLID, 
                         const double* x, Array<int>& facetLID)
  {
    if (cellDim==1)
      {
        static Array<int> facetLID(2);
        static Array<int> tmp(2);
        mesh.getFacetArray(cellDim, cellLID, 0, facetLID, tmp);
        Point A = mesh.nodePosition(facetLID[0]);
        Point B = mesh.nodePosition(facetLID[1]);
        if (A[0] < B[0]) return (x[0] >= A[0] && x[0] <= B[0]);
        else return (x[0] <= A[0] && x[0] >= B[0]);
      }
    else if (cellDim==2)
      {
        static Array<int> tmp(3);
        mesh.getFacetArray(cellDim, cellLID, 0, facetLID, tmp);
        const double* A = mesh.nodePositionView(facetLID[0]);
        const double* B = mesh.nodePositionView(facetLID[1]);
        const double* C = mesh.nodePositionView(facetLID[2]);
        /* first determine whether the three points of the triangle
         * are in ccw or cw order. */
        double sign = orient2D(A, B, C);
        if (sign > 0.0)
          {
            double s1 = orient2D(A, B, x);
            double s2 = orient2D(B, C, x);
            double s3 = orient2D(C, A, x);
            if (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) return true;
            return false;
          }
        else
          {
            double s1 = orient2D(A, C, x);
            double s2 = orient2D(C, B, x);
            double s3 = orient2D(B, A, x);
            if (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) return true;
            return false;
          }
      }
    else if (cellDim==3)
      {
        TEUCHOS_TEST_FOR_EXCEPT(true);
        return false; // -Wall
      }
    else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(cellDim<=0 || cellDim>3, std::runtime_error,
                           "invalid point dimension " << cellDim);
        return false; // -Wall
      }
  }

  double volume(const Mesh& mesh, int cellDim, int cellLID)
  {
    if (cellDim==1)
      {
        static Array<int> facetLID(2);
        static Array<int> tmp(2);
        mesh.getFacetArray(cellDim, cellLID, 0, facetLID, tmp);
        Point A = mesh.nodePosition(facetLID[0]);
        Point B = mesh.nodePosition(facetLID[1]);
        return fabs(A[0] - B[0]);
      }
    else if (cellDim==2)
      {
        static Array<int> facetLID(3);
        static Array<int> tmp(3);
        mesh.getFacetArray(cellDim, cellLID, 0, facetLID, tmp);
        const double* A = mesh.nodePositionView(facetLID[0]);
        const double* B = mesh.nodePositionView(facetLID[1]);
        const double* C = mesh.nodePositionView(facetLID[2]);
        return fabs(0.5*orient2D(A, B, C));
      }
    TEUCHOS_TEST_FOR_EXCEPT(true);
  }


  double orient2D(const double* A, const double* B, const double* x)
  {
    double acx, bcx, acy, bcy;

    acx = A[0] - x[0];
    bcx = B[0] - x[0];
    acy = A[1] - x[1];
    bcy = B[1] - x[1];
    return acx * bcy - acy * bcx;
  }

  int findEnclosingCell(const Mesh& mesh, int cellDim,
                        int initialGuessLID,
                        const double* x)
  {
    std::queue<int> Q;
    Set<int> repeats;
    static Array<int> facets;

    Q.push(initialGuessLID);

    

    while (!Q.empty())
      {
        int next = Q.front();
        Q.pop();
        if (repeats.contains(next)) continue;

        if (cellContainsPoint(mesh, cellDim, next, x, facets)) return next;
        repeats.put(next);
        
        std::list<int> neighbors;
        maximalNeighbors(mesh, cellDim, next, facets, neighbors);
        for (std::list<int>::const_iterator 
               i=neighbors.begin(); i!=neighbors.end(); i++)
          {
            Q.push(*i);
          }
      }
    return -1; // no containing cell found
  }

  Point pullback(const Mesh& mesh, int cellDim, int cellLID, const double* x)
  {
    int dim = cellDim;
    if (dim==2)
      {
        static Array<int> facetLID(3);
        static Array<int> tmp(3);
        mesh.getFacetArray(cellDim, cellLID, 0, facetLID, tmp);
        const double* A = mesh.nodePositionView(facetLID[0]);
        const double* B = mesh.nodePositionView(facetLID[1]);
        const double* C = mesh.nodePositionView(facetLID[2]);

        double bax = B[0] - A[0];
        double bay = B[1] - A[1];
        double cax = C[0] - A[0];
        double cay = C[1] - A[1];
        double delta = bax*cay - bay*cax;

        double xax = x[0] - A[0];
        double xay = x[1] - A[1];

        return Point( (cay*xax - cax*xay)/delta,
                      (-bay*xax + bax*xay)/delta );
          
      }
    else
      {
        TEUCHOS_TEST_FOR_EXCEPTION(cellDim != 2, std::runtime_error,
                           "invalid point dimension " << cellDim);
        return Point(); // -Wall
      }
  }

  void printCell(const Mesh& mesh, int cellLID)
  {
    if (mesh.spatialDim()==2)
      {
        static Array<int> facetLID(3);
        static Array<int> tmp(3);
        mesh.getFacetArray(mesh.spatialDim(), cellLID, 0, facetLID, tmp);
        Point A = mesh.nodePosition(facetLID[0]);
        Point B = mesh.nodePosition(facetLID[1]);
        Point C = mesh.nodePosition(facetLID[2]);
        std::cout << "{" << A << ", " << B << ", " << C << "}" << std::endl;
      }
  }

  void maximalNeighbors(const Mesh& mesh, int cellDim, int cellLID, 
                        const Array<int>& facetLID,
                        std::list<int>& rtn)
  {
    for (int f=0; f<facetLID.size(); f++)
      {
        Array<int> cofacets;
        mesh.getCofacets(0, facetLID[f], cellDim, cofacets);
        for (int c=0; c<cofacets.size(); c++)
          {
            if (cofacets[c] != cellLID) rtn.push_back(cofacets[c]);
          }
      }
  }


  int lookupEdgeLIDFromVerts(const Mesh& mesh, int v1, int v2)
  {
    Array<int> cofacetLIDs;
    mesh.getCofacets(0, v1, 1, cofacetLIDs);
    
    for (int c=0; c<cofacetLIDs.size(); c++)
      {
	int ori;
	for (int f=0; f<2; f++)
	  {
	    int v = mesh.facetLID(1, cofacetLIDs[c], 0, f, ori);
	    if (v == v2) return cofacetLIDs[c];
	  }
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
			       "edge (" << v1 << ", " << v2 << ") not found in mesh");
    return -1;
  }

  
  
}
