// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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
//
   
#ifndef BELOS_STATUS_TEST_MAXRESTARTS_HPP
#define BELOS_STATUS_TEST_MAXRESTARTS_HPP

/*!
  \file BelosStatusTestMaxRestarts.hpp
  \brief Belos::StatusTest class for specifying a maximum number of restarts.
*/

#include "BelosStatusTest.hpp"

/*! \class Belos::StatusTestMaxRestarts: 
    \brief A Belos::StatusTest class for specifying a maximum number of restarts.

    This implementation of the Belos::StatusTest base class tests the number of restarts performed
    against a maximum number allowed.  This status test is only valid for iterative methods which perform
    restarts, like GMRES.
*/

namespace Belos {

template <class ScalarType, class MV, class OP>
class StatusTestMaxRestarts: public StatusTest<ScalarType,MV,OP> {

 public:

   //! @name Constructor/Destructor.
  //@{ 

  //! Constructor
  StatusTestMaxRestarts(int maxIters);

  //! Destructor
  virtual ~StatusTestMaxRestarts() {};
  //@}

  //! @name Status methods
  //@{ 

  //! Check convergence status of the iterative solver: Unconverged, Converged, Failed.
  /*! This method checks to see if the convergence criteria are met using the current information from the 
    iterative solver.
  */
  StatusType CheckStatus(IterativeSolver<ScalarType,MV,OP> *iSolver );

  //! Return the result of the most recent CheckStatus call.
  StatusType GetStatus() const {return(status_);};

  //@}

  //! @name Reset methods
  //@{ 
  
  //! Resets the internal configuration to the initial state
  void Reset();

  //@}

  //! @name Accessor methods
  //@{ 

  //! Returns the maximum number of restarts set in the constructor.
  int GetMaxRestarts() const { return(maxRestarts_); };

  //! Returns the current number of restarts from the most recent StatusTest call.
  int GetNumRestarts() const { return(nRestarts_); };

  //@}

  //! @name Attribute methods
  //@{ 

  //! Indicates if residual vector is required by this convergence test: returns false for this class.
  bool ResidualVectorRequired() const { return(false); } ;
  //@}

  //! @name Print methods
  //@{ 

  //! Output formatted description of stopping test to output stream
  ostream& Print(ostream& os, int indent = 0) const;

  //@}
  

private:

  //! @name Private data members.
  //@{ 
  //! Maximum number of restarts allowed
  int maxRestarts_;

  //! Current number of restarts
  int nRestarts_;

  //! Status
  StatusType status_;
  //@}

};

  template <class ScalarType, class MV, class OP> 
  StatusTestMaxRestarts<ScalarType,MV,OP>::StatusTestMaxRestarts(int maxRestarts)
  {
    if (maxRestarts < 1)
      maxRestarts_ = 1;
    else
      maxRestarts_ = maxRestarts;
    
    nRestarts_ = 0;
    status_ = Undefined;
  }
  
  template <class ScalarType, class MV, class OP>
  StatusType StatusTestMaxRestarts<ScalarType,MV,OP>::CheckStatus(IterativeSolver<ScalarType,MV,OP> *iSolver )
  {
    status_ = Unconverged;
    nRestarts_ = iSolver->GetNumRestarts();
    if (nRestarts_ >= maxRestarts_)
      status_ = Failed;
    return status_;
  }
  
  template <class ScalarType, class MV, class OP>
  void StatusTestMaxRestarts<ScalarType,MV,OP>::Reset()
  {
    status_ = Undefined;
    nRestarts_ = 0;
  }

  template <class ScalarType, class MV, class OP>
  ostream& StatusTestMaxRestarts<ScalarType,MV,OP>::Print(ostream& os, int indent) const
  {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    this->PrintStatus(os, status_);
    os << "Number of Restarts = ";
    os << nRestarts_;
    os << ((nRestarts_ < maxRestarts_) ? " < " : ((nRestarts_ == maxRestarts_) ? " == " : " > "));
    os << maxRestarts_;
    os << endl;
    return os;
  }
  
} // end Belos namespace

#endif /* BELOS_STATUS_TEST_MAXRESTARTS_HPP */
