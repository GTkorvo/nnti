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


#ifndef TEUCHOS_XMLDEPENDENCYEXCEPTIONS_HPP_
#define TEUCHOS_XMLDEPENDENCYEXCEPTIONS_HPP_

/*! \file Teuchos_XMLDependencyExceptions.hpp
 * \brief A collection of Exceptions thrown
 * when converting Dependencys to and from
 * XML.
 */
#include <stdexcept>

namespace Teuchos {

/** \brief Thrown when no dependess of a dependency are specified
 * when converting the dependency from XML.
 */
class MissingDependeesException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingDependeesException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingDependeesException(const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when no dependents of a dependency are specified
 * when converting the dependency from XML.
 */
class MissingDependentsException : public std::logic_error{

public:

  /**
   * \brief Constructs an MissingDependentsException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  MissingDependentsException(const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when a Dependency has too many dependees specified in
 * its XML
 */
class TooManyDependeesException : public std::logic_error{

public:

  /**
   * \brief Constructs an TooManyDependeesException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  TooManyDependeesException(const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when a StringVisualDependency is being converted
 * from XML and no Values tag is found
 */
class ValuesTagMissingException : public std::logic_error{

public:

  /**
   * \brief Constructs an ValuesTagMissingException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  ValuesTagMissingException(const std::string& what_arg):std::logic_error(what_arg){}

};

/** \brief Thrown when a parent list can't be found in the
 * dependency sheet in which a dependency is trying to be
 * inserted.
 */
class ParentListNotFoundException : public std::logic_error{

public:

  /**
   * \brief Constructs an ParentListNotFoundException.
   *
   * @param what_arg The error message to be associated with this error.
   */
  ParentListNotFoundException(const std::string& what_arg):std::logic_error(what_arg){}

};


} // namespace Teuchos
#endif //TEUCHOS_XMLDEPENDENCYEXCEPTIONS_HPP_
