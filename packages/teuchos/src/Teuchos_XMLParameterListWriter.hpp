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

#ifndef Teuchos_XMLPARAMETERLISTWRITER_H
#define Teuchos_XMLPARAMETERLISTWRITER_H

/*! \file Teuchos_XMLParameterListWriter.hpp
    \brief Writes a ParameterList to an XML object
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_Utils.hpp"

namespace Teuchos
{

	/** \ingroup XML 
	 * \brief Writes a ParameterList to an XML object
	 */

	class TEUCHOS_LIB_DLL_EXPORT XMLParameterListWriter
		{
		public:
	  typedef std::map<const RCP<ParameterEntryValidator>, int> WriterValidatorIDMap;
	  typedef std::pair<const RCP<ParameterEntryValidator>, int> WriterValidatorIDPair;
      //! @name Constructors 
			//@{
      /** Construct a writer */
	  static const std::string validatorTagName;
      XMLParameterListWriter();
			//@}

      /** Write the given list to an XML object */
      XMLObject toXML(const ParameterList& p) const ;

	  static const std::string& getParameterListAspectsTagName(){
	    static const std::string parameterListAspectsTagName = "ParameterListAspects";
	    return parameterListAspectsTagName;
	  }

	  static const std::string& getParameterListsTagName(){
	    static const std::string parameterListsTagName = "ParameterLists";
	    return parameterListsTagName;
	  }

	  static const std::string& getParameterListTagName(){
	    static const std::string parameterListTagName = "ParameterList";
	    return parameterListTagName;
	  }

	  static const std::string& getNameAttributeName(){
	    static const std::string nameAttributeName = "name";
		return nameAttributeName;
	  }

	  static const std::string& getValidatorsTagName(){
	    static const std::string validatorsTagName = "Validators";
	    return validatorsTagName;
	  }

	  static const std::string& getValidatorIdAttributeName(){
	    static const std::string validatorIdAttributeName = "validatorid";
       return validatorIdAttributeName;
	  }

	
		private:

      /**
	   * \brief Write the given list to an XML object and record all the validators in it on a map.
	   */
      XMLObject convertParameterList(
	    const ParameterList& p,
		WriterValidatorIDMap& validatorIDMap,
		int& validatorIDCounter) const;
		};
}
#endif

