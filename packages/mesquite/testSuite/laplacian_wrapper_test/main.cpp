/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//   SUMMARY: 
//     USAGE:
//
// ORIG-DATE: 19-Feb-02 at 10:57:52
//  LAST-MOD: 14-Nov-02 at 17:33:19 by Thomas Leurent
//
//
// DESCRIPTION:
// ============
/*! \file main.cpp

describe main.cpp here

 */
// DESCRIP-END.
//

#ifdef USE_STD_INCLUDES
#include <iostream>
#else
#include <iostream.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MesquiteError.hpp"
#include "MeshSet.hpp"
#include "PlanarDomain.hpp"
// algorythms
#include "LaplacianIQ.hpp"


using namespace Mesquite;


#undef __FUNC__
#define __FUNC__ "main"
int main()
{
  Mesquite::MeshImpl *mesh = new Mesquite::MeshImpl;
  MsqError err;
  mesh->read_vtk("../../meshFiles/2D/VTK/square_quad_2.vtk", err);
     //create geometry: plane z=0, normal (0,0,1)
  Vector3D pnt(0,0,5);
  Vector3D s_norm(0,0,1);
  Mesquite::PlanarDomain msq_geom(s_norm, pnt, mesh);
    // initialises a MeshSet object
  MeshSet mesh_set1;
  mesh_set1.set_domain_constraint(&msq_geom, err); MSQ_CHKERR(err);
  mesh_set1.add_mesh(mesh, err); MSQ_CHKERR(err);
  
    // creates an intruction queue
  LaplacianIQ laplacian_smoother;
  
  mesh->write_vtk("original_mesh", err); MSQ_CHKERR(err);
  
    // launches optimization on mesh_set1
  laplacian_smoother.run_instructions(mesh_set1, err); MSQ_CHKERR(err);
  
  mesh->write_vtk("smoothed_mesh", err); MSQ_CHKERR(err);
}
