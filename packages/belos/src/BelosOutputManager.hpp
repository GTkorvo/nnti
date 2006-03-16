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

#ifndef BELOS_OUTPUT_MANAGER_HPP
#define BELOS_OUTPUT_MANAGER_HPP

/*!     \file BelosOutputManager.hpp
        \brief Class which manages the output and verbosity of the Belos solvers.
*/

#include "BelosConfigDefs.hpp"

/*!	\class Belos::OutputManager

	\brief Belos's basic output manager for sending information of select verbosity levels
	to the appropriate output stream.

	This output manager will remove the need for the solver or linear problem to know any information
	about the required output.  Calling <tt>doPrint( int vbLevel )</tt> will inform the solver if
	it is supposed to output the information corresponding to the verbosity level (\c vbLevel ).

	\author Michael Heroux and Heidi Thornquist
*/

namespace Belos {

  /*! \enum MsgType
    \brief Enumerated list of available message types recognized by the linear solvers.
  */
  
  enum MsgType 
    {
      Errors= 0,                  /*!< Errors [ always printed ] */
      Warnings = 0x1,             /*!< Internal warnings */
      IterationDetails = 0x2,     /*!< Approximate/exact residuals */
      OrthoDetails = 0x4,         /*!< Orthogonalization/orthonormalization details */
      FinalSummary = 0x8,         /*!< Final computational summary */
      TimingDetails = 0x10,       /*!< Timing details */
      Debug = 0x20                /*!< Debugging information */
    };
  
  template <class ScalarType>
  class OutputManager {
    
  public:
    
    //@{ \name Constructors/Destructor.
    
    //! Default constructor
    OutputManager();
    
    //! Basic constructor.
    OutputManager( int myID, int vbLevel = Belos::Errors, int printID = 0, ostream& os = std::cout );
    
    //! Destructor.
    virtual ~OutputManager() {};
    //@}
    
    //@{ \name Set methods.
    
    //! Set the output stream for this manager.
    void SetOStream( ostream& os ) { myOS_ = os; };
    
    //! Set the verbosity level for this manager.
    void SetVerbosity( int vbLevel ) { vbLevel_ = vbLevel; }; 
    
    //@}
    
    //@{ \name Get methods.
    
    //! Get the output stream for this manager.
    ostream& GetOStream() { return myOS_; };
    
    //@}
    
    //@{ \name Query methods.
    
    //! Find out whether we need to print out information for this message type.
    /*! This method is used by the solver to determine whether computations are
      necessary for this message type.
    */
    bool isVerbosity( MsgType type ) { return (( type == Belos::Errors ) || ( vbLevel_ & type )); }; 
    
    //! Find out whether this processor needs to print out information for this message type.
    /*! This method is used by the solver to determine whether this output stream has been
      selected to output the information for this message type.
    */
    bool isVerbosityAndPrint( MsgType type ) { return ( iPrint_ && isVerbosity( type )); }; 
    
    //! Find out whether information can be outputted through this output stream.
    bool doPrint( void ) { return (iPrint_); };
    
    //@}
    
  private:
    
    //@{ \name Undefined methods.
    
    //! Copy constructor.
    OutputManager( const OutputManager<ScalarType>& OM );
    
    //! Assignment operator.
    OutputManager<ScalarType>& operator=( const OutputManager<ScalarType>& OM );
    
    //@}
    
    int myID_, printID_;
    int vbLevel_;
    bool iPrint_;
    ostream& myOS_;	
  };
  
  template<class ScalarType>
  OutputManager<ScalarType>::OutputManager() :
    myID_(0),
    printID_(0),
    vbLevel_(0),
    iPrint_(true),
    myOS_(std::cout)
  {
  }
  
  template<class ScalarType>
  OutputManager<ScalarType>::OutputManager( int myID, int vbLevel, int printID, ostream& os ) :
    myID_(myID),
    printID_(printID),
    vbLevel_(vbLevel),
    iPrint_(myID == printID),
    myOS_(os)
  {
  }
  
} // end Belos namespace

#endif

// end of file BelosOutputManager.hpp
