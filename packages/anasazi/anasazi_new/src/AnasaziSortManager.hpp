// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_SORT_MANAGER_HPP
#define ANASAZI_SORT_MANAGER_HPP

/*!     \file AnasaziSortManager.hpp
        \brief Class which manages the sorting of approximate eigenvalues computed by Anasazi solvers.
*/

/*!    \class Anasazi::SortManager
       \brief Anasazi's templated pure virtual class for managing the sorting of 
       approximate eigenvalues computed by the eigensolver.

       A concrete implementation of this class is necessary.  The user can create
       their own implementation if those supplied are not suitable for their needs.

       \author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziReturnType.hpp"
#include "AnasaziEigensolver.hpp"

namespace Anasazi {

  template<class STYPE, class MV, class OP>
  class SortManager {
    
  public:
    
    //! Default Constructor
    SortManager() {};

    //! Destructor
    virtual ~SortManager() {};

    //! Sort the vector of eigenvalues, optionally returning the permutation vector.
    /**
       @param solver [in] Eigensolver that is calling the sorting routine

       @param n [in] Size of the array

       @param evals [in/out] Array of length n containing the eigenvalues to be sorted

       @param perm [out] Vector of length n to store the permutation (optional)

       @return Returns the status of the sorting routine [ Undefined by default ] 
    */
    virtual ReturnType sort(Eigensolver<STYPE,MV,OP>* solver, int n, STYPE *evals, std::vector<int> *perm = 0) const { return Undefined; };
    
    //! Sort the vectors of eigenpairs, optionally returning the permutation vector.
    /**
       @param solver [in] Eigensolver that is calling the sorting routine

       @param n [in] Size of the array

       @param r_evals [in/out] Array of length n containing the real part of the eigenvalues to be sorted 

       @param i_evals [in/out] Array of length n containing the imaginary part of the eigenvalues to be sorted 

       @param perm [out] Vector of length n to store the permutation (optional)

       @return Returns the status of the sorting routine [ Undefined by default ] 
    */
    virtual ReturnType sort(Eigensolver<STYPE,MV,OP>* solver, int n, STYPE *r_evals, STYPE *i_evals, std::vector<int> *perm = 0) const { return Undefined; };
    
  };
  
}

#endif // ANASAZI_SORT_MANAGER_HPP

