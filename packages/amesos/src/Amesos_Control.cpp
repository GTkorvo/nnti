// @HEADER
// ***********************************************************************
//
//                Amesos: Direct Sparse Solver Package
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

#include "Amesos_Control.h"
void Amesos_Control::SetControlParameters(Teuchos::ParameterList &ParameterList) {

  // add zero to diagonal if diagonal element is not present
  // - some solvers choke on matrices with zero elements on the diagonal
  if( ParameterList.isParameter("AddZeroToDiag") )
    AddZeroToDiag_ = ParameterList.get("AddZeroToDiag", AddZeroToDiag_);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);

  // threshold for determining if refactorize worked OK
  // UNUSED at present - KSS June 2004
  if( ParameterList.isParameter("RcondThreshold") )
    rcond_threshold_ = ParameterList.get("RcondThreshold", 1e-12);

  // define how many processes to use in the ScaLAPACK factor and solve
  // if (-1), a heuristic is used to determine the number of processes to use 
  if( ParameterList.isParameter("MaxProcs") )
    MaxProcesses_ = ParameterList.get("MaxProcs",MaxProcesses_);
  
  // Matrix property, defined internally in Amesos_Mumps as an integer,
  // whose value can be:
  // - 0 : general unsymmetric matrix;
  // - 1 : SPD;
  // - 2 : general symmetric matrix.
  if( ParameterList.isParameter("MatrixProperty") ) {
    string MatrixProperty;
    MatrixProperty = ParameterList.get("MatrixProperty",MatrixProperty);
    if( MatrixProperty == "SPD" )
      MatrixProperty_ = 1;
    else if( MatrixProperty == "symmetric" ) 
      MatrixProperty_ = 2;
    else if( MatrixProperty == "general" )   
      MatrixProperty_ = 0;
    else {
      //      AMESOS_CHK_ERR( -1 ) ;
      //      if ( verbose_ ) cerr << "Amesos : ERROR" << endl 
      //	     << "Amesos : MatrixProperty value not recognized ("
      //	     << MatrixProperty << ")" << endl;
    }
  }

  // scaling method: 0: none, 1: use method's default, 2: use
  // the method's 1st alternative, 3: etc.
  if( ParameterList.isParameter("ScaleMethod") )
    ScaleMethod_ = ParameterList.get("ScaleMethod", ScaleMethod_);


  if( ParameterList.isParameter("Reindex") )
    Reindex_ = ParameterList.get("Reindex", Reindex_);

}
