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

#ifndef ANASAZI_OUTPUT_MANAGER_HPP
#define ANASAZI_OUTPUT_MANAGER_HPP

/*!     \file AnasaziOutputManager.hpp
        \brief Class which manages the output and verbosity of the Anasazi solvers.
*/

#include "AnasaziConfigDefs.hpp"

/*!	\class Anasazi::OutputManager

	\brief Anasazi's basic output manager for sending information of select verbosity levels
	to the appropriate output stream.

	This output manager will remove the need for the solver or linear problem to know any information
	about the required output.  Calling <tt>doOutput( int vbLevel )</tt> will inform the solver if
	it is supposed to output the information corresponding to the verbosity level (\c vbLevel ).

	\author Michael Heroux and Heidi Thornquist
*/

namespace Anasazi {

template <class ScalarType>
class OutputManager {

	public:

	//@{ \name Constructors/Destructor.

	//! Default constructor
	OutputManager();

	//! Copy constructor.
	OutputManager( const OutputManager<ScalarType>& OM );

	//! Basic constructor.
	OutputManager( int myID, int vbLevel = 0, int printID = 0, ostream& os = cout );

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

	//! Find out whether we need to print out information for this verbosity level.
	bool doOutput( int vbLevel ) { return (( IPrint_ && vbLevel_ > vbLevel ) ? true : false ); }; 

	//@}

	private:

	int myID_, printID_;
	int vbLevel_;
	bool IPrint_;
	ostream& myOS_;	
};

template<class ScalarType>
OutputManager<ScalarType>::OutputManager() :
	myID_(0),
	printID_(0),
	vbLevel_(0),
	IPrint_(true),
	myOS_(cout)
{
}

template<class ScalarType>
OutputManager<ScalarType>::OutputManager( const OutputManager<ScalarType>& OM ) :
	myID_(OM.myID_),
	printID_(OM.printID_),
	vbLevel_(OM.vbLevel_),
	IPrint_(OM.IPrint_),
	myOS_(OM.myOS_)
{
}	 

template<class ScalarType>
OutputManager<ScalarType>::OutputManager( int myID, int vbLevel, int printID, ostream& os ) :
	myID_(myID),
	printID_(printID),
	vbLevel_(vbLevel),
	IPrint_(myID == printID),
	myOS_(os)
{
}

} // end Anasazi namespace

#endif

// end of file AnasaziOutputManager.hpp
