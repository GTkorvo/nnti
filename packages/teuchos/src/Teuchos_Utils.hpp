// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_UTILS_H
#define TEUCHOS_UTILS_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos
{
  using std::string;

  /**\ingroup Utilities
   *
   * Numerical constants, etc.
   */
  class Utils
    {
    public:
      /**
       * print a description of the current build
       */
      static void aboutBuild();

      /**
       * Set a number to zero if it is less than chopVal.
       */
      static double chop(const double& x);

      /**
       * write a real as a string
       */
      static string toString(const double& x);

      /**
       * write an int as a string
       */
      static string toString(const int& x);

      /**
       * IEEE positive infinity
       */
      static double infinity() {return HUGE_VAL;}

      /**
       * IEEE negative infinity
       */
      static double negativeInfinity() {return -HUGE_VAL;}

      /**
       * pi.
       */
      static double pi() {return M_PI;}

      /**
       * Get the chopping value, below which numbers are considered to be zero
       */
      static double getChopVal() {return chopVal_;}
      /**
       * Set the chopping value, below which numbers are considered to be zero
       */
      static void setChopVal(double chopVal) {chopVal_ = chopVal;}

    private:
      static double chopVal_;
    };

  /** \relates Utils */
  inline string toString(const int& x) {return Utils::toString(x);}

  /** \relates Utils */
  inline string toString(const double& x) {return Utils::toString(x);}

  /** \relates Utils */
  inline string toString(const string& x) {return x;}

}

#endif


