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

#ifndef SUNDANCE_WATCHFLAG_H
#define SUNDANCE_WATCHFLAG_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_XMLObject.hpp"
#include <string>

namespace SundanceUtils
{
  /** 
   * WatchFlag
   */
  class WatchFlag
    {
    public:
      /** */
      WatchFlag(const std::string& name="");

      /** */
      const std::string& name() const {return name_;}

      /** */
      void activate() ;

      /** */
      void deactivate() ;

      /** */
      bool isActive() const ;

      /** */
      bool operator<(const WatchFlag& other) const
        {return name() < other.name();}

      /** */
      XMLObject toXML() const ;

      /** */
      int evalVerb() const {return evalVerbosity_;}

      /** */
      int setupVerb() const {return setupVerbosity_;}

      /** */
      void setEvalVerb(int v) {evalVerbosity_=v;}

      /** */
      int setSetupVerb(int v) {setupVerbosity_=v;}

    private:
      std::string name_;

      int evalVerbosity_;

      int setupVerbosity_;

      static Map<std::string, bool>& isActiveMap()
        {
          static Map<std::string, bool> rtn;
          return rtn;
        }

    };
}

#endif
