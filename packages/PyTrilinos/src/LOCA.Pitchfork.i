// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%module(package="PyTrilinos.LOCA") Pitchfork

%{
#include "LOCA_Extended_MultiAbstractGroup.H"
#include "LOCA_BorderedSystem_AbstractGroup.H"
#include "LOCA_MultiContinuation_ExtendedGroup.H"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_AbstractStrategy.H"

#include "LOCA_MultiContinuation_FiniteDifferenceGroup.H"
#include "LOCA_TurningPoint_MooreSpence_FiniteDifferenceGroup.H"
#include "LOCA_TurningPoint_MinimallyAugmented_FiniteDifferenceGroup.H"

#include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"

// Local includes
#define NO_IMPORT_ARRAY
#include "numpy_include.h"

// Namespace flattening
using Teuchos::RCP;
%}

// Standard exception handling
%include "exception.i"

// Ignore/renames
%ignore *::operator=;

// Trilinos module imports
%import "Teuchos.i"

%import "LOCA.TurningPoint.i"

%rename(MooreSpence_AbstractGroup) LOCA::Pitchfork::MooreSpence::AbstractGroup;
%rename(MinimallyAugmented_AbstractGroup) LOCA::Pitchfork::MinimallyAugmented::AbstractGroup;

%include "LOCA_Pitchfork_MooreSpence_AbstractGroup.H"
%include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
