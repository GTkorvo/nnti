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

#ifndef SUNDANCE_MULTISET_H
#define SUNDANCE_MULTISET_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"
#include <set>

#ifndef DOXYGEN_DEVELOPER_ONLY


namespace SundanceUtils
{
  using namespace Teuchos;

  /** */
  template<class Key>
    class MultiSet : public std::multiset<Key>
    {
    public:
      /** */
      MultiSet() : std::multiset<Key>() {;}

      /** */
      bool contains(const Key& key) const {return find(key) != end();}

      /** */
      void put(const Key& key) {insert(key);}

      /** */
      Array<Key> elements() const ;

      /** */
      ostream& toStream(ostream& os) const ;

      /** */
      MultiSet<Key> merge(const MultiSet<Key>& other) const ;

      /** */
      string toString() const ;
    };


  template<class Key> inline
    Array<Key> MultiSet<Key>::elements() const
    {
      Array<Key> rtn;

      typename MultiSet<Key>::const_iterator iter;

      for (iter=begin(); iter != end(); iter++)
        {
          rtn.append(*iter);
        }
      return rtn;
    }

  template<class Key> inline
  MultiSet<Key> MultiSet<Key>::merge(const MultiSet<Key>& other) const
    {
      MultiSet<Key> rtn = *this;

      typename MultiSet<Key>::const_iterator iter;

      for (iter=other.begin(); iter != other.end(); iter++)
        {
          rtn.put(*iter);
        }
      return rtn;
    }

  template<class Key> inline
    ostream& MultiSet<Key>::toStream(ostream& os) const
    {
      typename MultiSet<Key>::const_iterator iter;

      unsigned int k = 0;
      os << "{";
      for (iter=begin(); iter != end(); iter++, k++)
        {
          os << *iter;
          if (k<(size()-1)) os << ", ";
        }
      os << "}";

      return os;
    }

  template<class Key> inline
    string MultiSet<Key>::toString() const
    {
      ostringstream os;
      os << *this;
      return os.str();
    }

}

namespace std
{
/** \relates MultiSet */
  template<class Key> inline
    ostream& operator<<(ostream& os, const SundanceUtils::MultiSet<Key>& m)
    {return m.toStream(os);}
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif

