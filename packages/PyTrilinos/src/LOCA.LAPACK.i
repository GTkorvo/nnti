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

%module(package="PyTrilinos.LOCA") LAPACK

%{
// LOCA includes
#include "LOCA.H"
#include "LOCA_LAPACK.H"
%}

// Ignore/renames
%rename(Print) *::print() const;
%ignore *::operator=;
%ignore operator<<(ostream& stream, const NOX::LAPACK::Vector& v);

// Import LOCA interface
%import "LOCA_Abstract.i"

// Import NOX_LAPACK headers
%import "NOX_LAPACK_Interface.H"
%import "NOX_LAPACK_Group.H"

// LOCA interface includes
%include "LOCA_LAPACK_Interface.H"
%include "LOCA_LAPACK_Group.H"
