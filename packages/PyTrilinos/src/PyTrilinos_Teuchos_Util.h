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

#ifndef PYTRILINOS_TEUCHOS_UTIL_H
#define PYTRILINOS_TEUCHOS_UTIL_H

// Include python headers
#include "Python.h"
#ifdef HAVE_INTTYPES_H
#undef HAVE_INTTYPES_H
#endif

// Include Teuchos::ParameterList prototypes
#include "Teuchos_ParameterList.hpp"

// ****************************************************************** //

// The Teuchos::ParameterList class can support parameters of any
// type, but the python wrappers need a subset of concrete parameter
// types supported a-priori.  Since this subset shows up in many
// places, the following functions are provided to handle conversion
// between PyObjects and the supported C / C++ types in a single
// place.  The following type conversions are supported:

//    Python                           C / C++
//    ---------------------      -------------
//    bool       	    <--> bool
//    int        	    <--> int
//    float      	    <--> double
//    string     	    <--> std::string
//    string     	    <--  char *
//    dict       	     --> ParameterList
//    wrapped ParameterList <--> ParameterList

//    Note: The python None type is unsupported and used for error
//    reporting in getPythonParameter().

// To convert the wrapped ParameterList class, certain SWIG functions
// are used that only exist within the wrapper code generated by swig.
// Therefore, this code only compiles correctly within the context of
// a swig interface file.

// ****************************************************************** //

namespace PyTrilinos
{

// ****************************************************************** //

// The following enumeration is used as an optional flag in certain
// routines to indicate how the routine is supposed to react to
// illegal parameters.

enum ResponseToIllegalParameters {raiseError,
				  ignore,
				  storeNames };

// ****************************************************************** //

// The setPythonParameter() function takes a Teuchos::ParameterList, a
// string name and a python object as input.  An attempt is made to
// convert the python object to a supported C / C++ type.  If
// successful, the Teuchos::ParameterList set() method is called using
// the given name and converted object, and true is returned.  If
// unsuccessful (ie, the python object represents an unsupported
// type), false is returned.  Typically, a python error will NOT be
// set except in the case when the python object is a dictionary with
// unsupported values.

bool setPythonParameter(Teuchos::ParameterList & plist,
			const std::string      & name,
			PyObject               * value);

// **************************************************************** //

// The getPythonParameter() function is the get counterpart to
// setPythonParameter().  It takes a Teuchos::ParameterList and string
// name as input.  If the requested parameter name does not exist, a
// new reference to None is returned (a type that is guaranteed not to
// be supported by setPythonParameter()).  If the name exists and its
// type is supported, it is returned as a new reference to a python
// object.  If the name exists, but the type is not supported, NULL is
// returned, to indicate an error.  All returned python object
// pointers are new references.  This function is coded in such a way
// that the Teuchos::ParameterList "used" flags are not altered.

PyObject * getPythonParameter(const Teuchos::ParameterList & plist,
			      const std::string            & name);

// **************************************************************** //

// Function isEquivalent() is a utility for determining whether a
// python dictionary is functionally equivalent to a
// Teuchos::ParameterList.  It supports interpreting
// Teuchos::ParameterList sublists as nested python dictionaries, so
// it calls itself recursively.

bool isEquivalent(PyObject                     * dict,
		  const Teuchos::ParameterList & plist);

// **************************************************************** //

// Function updatePyDictWithParameterList() takes all of the entries
// in a Teuchos::ParameterList and updates the given python dictionary
// to reflect the same values.  If the given python object is not a
// dictionary, or any of the Teuchos::ParameterList entries are of
// unsupported type, the function returns false.

bool updatePyDictWithParameterList(PyObject                     * dict,
				   const Teuchos::ParameterList & plist,
				   ResponseToIllegalParameters    flag=raiseError);

// **************************************************************** //

// Function updateParameterListWithPyDict() takes all of the entries
// in a python dictionary and updates the given Teuchos::ParameterList
// to reflect the same values.  If the given python object is not a
// dictionary, or if any of the dictionary keys are not strings, or if
// any of the dictionary values are of unsupported type, then the
// function returns false.

bool updateParameterListWithPyDict(PyObject                  * dict,
				   Teuchos::ParameterList    & plist,
				   ResponseToIllegalParameters flag=raiseError);

// **************************************************************** //

// Function synchronizeParameters() is a function for bringing the
// given python dictionary and Teuchos::ParameterList into synch with
// each other.  If a parameter exists for both the
// Teuchos::ParameterList and the python dictionary, the
// Teuchos::ParameterList takes precedence.  If the function returns
// false, it means the given PyObject was not a dictionary or the
// Teuchos::ParameterList or python dictionary had at least one value
// of an unsupported type.

bool synchronizeParameters(PyObject                  * dict,
			   Teuchos::ParameterList    & plist,
			   ResponseToIllegalParameters flag=raiseError);

// **************************************************************** //

// Function pyDictToNewParameterList is a helper function that takes a
// python dictionary and returns a pointer to an equivalent, new
// Teuchos::ParameterList.  If dict is not a python dictionary, or dict is
// not a valid dictionary (non-string keys or unsupported value
// types) then the function returns NULL.

Teuchos::ParameterList * pyDictToNewParameterList(PyObject                  * dict,
						  ResponseToIllegalParameters flag=raiseError);

// **************************************************************** //

// Function parameterListToNewPyDict is a helper function that takes a
// Teuchos::ParameterList and returns a pointer to an equivalent, new
// python dictionary.  If the Teuchos::ParameterList contains entries
// of invalid type, then a python error is raised and NULL is
// returned.

PyObject * parameterListToNewPyDict(const Teuchos::ParameterList & plist,
				    ResponseToIllegalParameters    flag=raiseError);

}    // Namespace PyTrilinos

#endif // PYTRILINOS_TEUCHOS_UTIL_H
