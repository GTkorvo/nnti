/*
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
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"
#include "Teuchos_ParameterEntryXMLConverterDB.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"
#include <iostream>



/**
* When parameterlists get written or read from XML, any validators 
* associated with a parameter entry go with it. If you're using a
* validator that isn't part of the standard set of validators, then
* you'll get an error. You need to add a converter for it.
*/
int main(int argc, char* argv[])
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  Teuchos::ParameterList My_List;

  //Basic data types
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<short> >
    solverValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<short>(
        Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )
        ,"Solver"
        )
      );
  My_List.set(
    "Solver"
    ,"GMRES" 
    ,"The type of solver to use."
    ,solverValidator
    );

  //By default, there is no validator converter for a StringToIntegralParameterEntryValidator
  //templated on type short. So lets add one with the following convience macro found
  //in ValidatorXMLConverterDB.hpp
  TEUCHOS_ADD_STRINGTOINTEGRALCONVERTER(short)

  //Of course if you have a completly custom validator you'll need to do a little more
  //You'll have to make you're own converter for it by subclassing the ValidatorXMLConverter
  //class and then use the TEUCHOS_ADD_VALIDATOR_CONVERTER convience macro.

  
   //Now we'll write it out to xml.
  Teuchos::writeParameterListToXmlFile(My_List, "My_List.xml");
  //Then read it in to a new list.
 
  Teuchos::writeParameterListToXmlOStream(My_List, std::cout);
  Teuchos::RCP<Teuchos::ParameterList> readIn = Teuchos::getParametersFromXmlFile("My_List.xml");

  std::cout << *readIn;

  /**
   * Final Notes:
   * StandardTemplatedParameterConverter should suit most your needs. Buf if for some reason you
   * don't feel like overrideing the inseration and extraction operators, you can allways subclass
   * the ParameterEntryXMLConverter class and do your own thing.
   */
  return 0;
}
