// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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

#ifndef BJBRIGHTOPERATORSOURCE_H
#define BJBRIGHTOPERATORSOURCE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{
  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * BJBRightOperatorSource
   * X = C - B D^{-1} B^T, where D = diag(conv-diff op)
   *
   */

  // Need a better name for this class.

  class BJBRightOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      // BJBRightOperatorSource(TSFLinearOperator& S);

      /** ctor (we make Ap = C from saddle matrix S) */
      BJBRightOperatorSource(TSFLinearOperator& S);

      BJBRightOperatorSource(TSFLinearOperator& S,
                             TSFReal alp, TSFReal beta);
      /** ctor (we are given Ap) */
      BJBRightOperatorSource(TSFLinearOperator& S,
                             TSFLinearOperator& Ap);

      BJBRightOperatorSource(TSFLinearOperator& S,
                             TSFLinearOperator& Ap,
			     TSFReal alp, TSFReal beta);
      //get beta, alp routine!
      /** virtual destructor */
      virtual ~BJBRightOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getAp() const;

      TSFLinearOperator getDinv() const;

      TSFReal getBeta() const;

      TSFReal getAlpha() const;

      bool getScal() const;

      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      mutable TSFReal alp_;
      mutable TSFReal beta_;
      mutable TSFLinearOperator Dinv_;
      mutable TSFLinearOperator Ap_;
      mutable bool hasAp_;
      mutable bool hasscal_;
      };


}


#endif
