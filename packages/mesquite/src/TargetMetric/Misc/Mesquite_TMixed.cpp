/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2009 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA

    (2009) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file TMixed.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "Mesquite_TMetricBarrier.hpp"
#include "Mesquite_TMixed.hpp"
#include "Mesquite_MsqMatrix.hpp"
#include "Mesquite_MsqError.hpp"
#include <sstream>

namespace MESQUITE_NS {

std::string TMixed::get_name() const
{
  std::ostringstream name;
  name << "2D:" << mu2D->get_name() << "; "
       << "3D:" << mu3D->get_name();
  return name.str();
}

TMixed::~TMixed() {}

bool TMixed::evaluate( const MsqMatrix<2,2>& T, 
                       double& result, 
                       MsqError& err )
{
  bool rval = mu2D->evaluate( T, result, err );
  MSQ_ERRZERO(err);
  return rval;
}

bool TMixed::evaluate( const MsqMatrix<3,3>& T, 
                       double& result, 
                       MsqError& err )
{
  bool rval = mu3D->evaluate( T, result, err );
  MSQ_ERRZERO(err);
  return rval;
}

bool TMixed::evaluate_with_grad( const MsqMatrix<2,2>& T,
                                 double& result,
                                 MsqMatrix<2,2>& deriv_wrt_T,
                                 MsqError& err )
{
  bool rval = mu2D->evaluate_with_grad( T, result, deriv_wrt_T, err );
  MSQ_ERRZERO(err);
  return rval;
}

bool TMixed::evaluate_with_grad( const MsqMatrix<3,3>& T,
                                 double& result,
                                 MsqMatrix<3,3>& deriv_wrt_T,
                                 MsqError& err )
{
  bool rval = mu3D->evaluate_with_grad( T, result, deriv_wrt_T, err );
  MSQ_ERRZERO(err);
  return rval;
}

bool TMixed::evaluate_with_hess( const MsqMatrix<2,2>& T,
                                 double& result,
                                 MsqMatrix<2,2>& deriv_wrt_T,
                                 MsqMatrix<2,2> second_wrt_T[3],
                                 MsqError& err )
{
  bool rval = mu2D->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );
  MSQ_ERRZERO(err);
  return rval;
}

bool TMixed::evaluate_with_hess( const MsqMatrix<3,3>& T,
                                 double& result,
                                 MsqMatrix<3,3>& deriv_wrt_T,
                                 MsqMatrix<3,3> second_wrt_T[3],
                                 MsqError& err )
{
  bool rval = mu3D->evaluate_with_hess( T, result, deriv_wrt_T, second_wrt_T, err );
  MSQ_ERRZERO(err);
  return rval;
}


} // namespace MESQUITE_NS
