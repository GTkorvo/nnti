// @HEADER 
// ***********************************************************************
// 
//         Optika: A Tool For Developing Parameter Obtaining GUIs
//                Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, with Sandia Corporation, the 
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Kurtis Nusbaum (klnusbaum@gmail.com) 
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_XMLParameterListExceptions.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
#include "Teuchos_ValidatorXMLConverterDB.hpp"
#include "Teuchos_StandardValidatorXMLConverters.hpp"


namespace Teuchos {

class UNDEFINED_PARAMETERENTRY_VALIDATOR : public ParameterEntryValidator
{
  
  public:

  void printDoc(const std::string& docString, std::ostream& out) const {}

  ValidStringsList validStringValues() const{
    return rcp(new Array<std::string>(1,""));
  }
  
  void validate(
    ParameterEntry  const& entry,
    std::string const& paramName,
    std::string const& sublistName
    ) const {}

  const std::string getXMLTypeName() const{
    return "UNDEFINEDTYPE";
  }

};

TEUCHOS_UNIT_TEST(Teuchos_Validator, exceptionTests)
{
  ValidatorXMLConverterDB::printKnownConverters(out);
  out << std::endl;

  UNDEFINED_PARAMETERENTRY_VALIDATOR badValidator;
  TEST_THROW(ValidatorXMLConverterDB::getConverter(badValidator), CantFindValidatorConverterException);

  TEST_THROW(RCP<ParameterList>
    missingValidatorList = getParametersFromXmlFile("MissingValidator.xml"),
    MissingValidatorDefinitionException);
 
  TEST_THROW(RCP<ParameterList>
    missingPrototypeList = getParametersFromXmlFile("MissingPrototypeValidator.xml"),
	MissingValidatorDefinitionException);

  TEST_THROW(RCP<ParameterList>
    conflicitingValiIdsList = getParametersFromXmlFile("ConflictingValidatorIDs.xml"),
    DuplicateValidatorIDsException);

  TEST_THROW(RCP<ParameterList>
    stringValidatorBadTagList = getParametersFromXmlFile("StringValidatorBadTag.xml"),
    BadTagException);

  TEST_THROW(RCP<ParameterList>
    stringValidatorBadTagList = getParametersFromXmlFile("StringToIntegralValidatorBadTag.xml"),
    BadTagException);

  #ifdef HAVE_TEUCHOS_DEBUG

  StringValidatorXMLConverter stringConverter;
  AnyNumberValidatorXMLConverter anyNumberConverter;
  XMLParameterListWriter::ValidatorIDsMap writerDummyMap;
  XMLParameterListReader::ValidatorIDsMap readerDummyMap;
  RCP<AnyNumberParameterEntryValidator> anyNumberValidator = 
    anyNumberParameterEntryValidator();
  writerDummyMap.insert(XMLParameterListWriter::ValidatorIDsMap::value_type(
    anyNumberValidator, 0));
  TEST_THROW(
    stringConverter.fromValidatortoXML(anyNumberValidator, writerDummyMap), 
    BadValidatorXMLConverterException);
  XMLObject anyNumberXML = 
    anyNumberConverter.fromValidatortoXML(anyNumberValidator, writerDummyMap);
  TEST_THROW(
    stringConverter.fromXMLtoValidator(anyNumberXML, readerDummyMap), 
    BadValidatorXMLConverterException);

  #endif

}

TEUCHOS_UNIT_TEST(Teuchos_Validator, fileNameValidatorConverter)
{
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "non default";

  RCP<FileNameValidator> defaultValidator =
    rcp(new FileNameValidator);
  RCP<FileNameValidator> nonDefaultValidator =
    rcp(new FileNameValidator(true));
  ParameterList myList("FileName Validator List");
  myList.set("default", "", "parameter for default validator",
    defaultValidator);
  myList.set("non default", "blah.txt", "parameter for non default validator",
    nonDefaultValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const FileNameValidator> readinDefault =
    rcp_dynamic_cast<const FileNameValidator>(
      readInPL->getEntry(defaultParameterName).validator(), true);
  TEST_EQUALITY(readinDefault->fileMustExist(), defaultValidator->fileMustExist());

  RCP<const FileNameValidator> readinNonDefault =
    rcp_dynamic_cast<const FileNameValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator(), true);
  TEST_EQUALITY(readinNonDefault->fileMustExist(), nonDefaultValidator->fileMustExist());
}


TEUCHOS_UNIT_TEST(Teuchos_Validator, stringValidatorConverter)
{
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "non default";

  RCP<StringValidator> nonDefaultValidator = rcp(
    new StringValidator(tuple<std::string>("value1", "cheese", "kurtis", "is", "awesome")));
  ParameterList myList("String Validator List");
  myList.set("non default", "kurtis", "parameter for non default validator",
    nonDefaultValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const StringValidator> readinNonDefault =
    rcp_dynamic_cast<const StringValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator(), true);
  TEST_COMPARE_ARRAYS(*(readinNonDefault->validStringValues()),
    *(nonDefaultValidator->validStringValues()));
}


TEUCHOS_UNIT_TEST(Teuchos_Validator, anynumberValidatorConverter)
{
  std::string xmlFileName = "AnyNumberValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string nonDefaultParameterName = "preferred and accepted";
  RCP<AnyNumberParameterEntryValidator> defaultValidator =
    rcp(new AnyNumberParameterEntryValidator());
  AnyNumberParameterEntryValidator::AcceptedTypes acceptedTypes;
  acceptedTypes.allowDouble(false);
  RCP<AnyNumberParameterEntryValidator> nonDefaultValidator =
    rcp(
      new AnyNumberParameterEntryValidator(
        AnyNumberParameterEntryValidator::PREFER_INT,
        acceptedTypes
        )
      );

  ParameterList myList("AnyNumberValidatorList");
  myList.set(defaultParameterName, 10.0,
    "A parameter with the default AnyNumberValidator on it", defaultValidator);
  myList.set(nonDefaultParameterName, 1, 
    "A prameter with an AnyNumberValidator on it that has the preferred and accepted types differnet from the default",
    nonDefaultValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);
  
  RCP<const AnyNumberParameterEntryValidator> readinDefaultValidator =
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      readInPL->getEntry(defaultParameterName).validator(), true);
  TEST_EQUALITY(readinDefaultValidator->isDoubleAllowed(),
    defaultValidator->isDoubleAllowed());
  TEST_EQUALITY(readinDefaultValidator->isIntAllowed(),
    defaultValidator->isIntAllowed());
  TEST_EQUALITY(readinDefaultValidator->isStringAllowed(),
    defaultValidator->isStringAllowed());
  TEST_EQUALITY(readinDefaultValidator->getPreferredType(),
    defaultValidator->getPreferredType());

  RCP<const AnyNumberParameterEntryValidator> readinNonDefaultValidator =
    rcp_dynamic_cast<const AnyNumberParameterEntryValidator>(
      readInPL->getEntry(nonDefaultParameterName).validator(), true);
  TEST_EQUALITY(readinNonDefaultValidator->isDoubleAllowed(),
    nonDefaultValidator->isDoubleAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->isIntAllowed(),
    nonDefaultValidator->isIntAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->isStringAllowed(),
    nonDefaultValidator->isStringAllowed());
  TEST_EQUALITY(readinNonDefaultValidator->getPreferredType(),
    nonDefaultValidator->getPreferredType());
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, EnhancedNumberValidatorConverter, T)
{
  std::string xmlFileName = TypeNameTraits<T>::name() + "EnhancedValidatorList.xml";
  std::string defaultParameterName = "default";
  std::string minmaxParameterName = "min max";
  ParameterList myList;
  RCP<EnhancedNumberValidator< T > > defaultValidator =
    rcp( new EnhancedNumberValidator< T >());
  RCP<EnhancedNumberValidator< T > > minMaxValidator =
    rcp( new EnhancedNumberValidator< T >(0,10));
  myList.set(defaultParameterName, ( T )6, "parameter with default validator",
    defaultValidator);
  myList.set(minmaxParameterName, ( T )10, "parameter with min and max validator",
    minMaxValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->getMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->getMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->getMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->getMax()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->getStep()
    ,
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->getStep()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(
        defaultParameterName).validator(), true)->getPrecision(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(
        defaultParameterName).validator(), true)->getPrecision()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->hasMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->hasMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(defaultParameterName).validator(), true)->hasMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(defaultParameterName).validator(), true)->hasMax()
  );

  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->getMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->getMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->getMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->getMax()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->getStep(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->getStep()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(
        minmaxParameterName).validator(), true)->getPrecision(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(
        minmaxParameterName).validator(), true)->getPrecision()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->hasMin(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->hasMin()
  );
  TEST_EQUALITY(
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      readInPL->getEntry(minmaxParameterName).validator(), true)->hasMax(),
    rcp_dynamic_cast<const EnhancedNumberValidator< T > >(
      myList.getEntry(minmaxParameterName).validator(), true)->hasMax()
  );

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, NumberArrayValidatorConverterTest, T)
{
  std::string arrayParameterName = "array";
  ParameterList myList;

  const T arrayValidatorLen = as<T>(11);
  RCP<ArrayNumberValidator< T > > arrayValidator =
    rcp(new ArrayNumberValidator< T >(
      rcp(new EnhancedNumberValidator<T>(as<T>(0), arrayValidatorLen))));
  myList.set(arrayParameterName,
    Array< T >(4, 10), "array parameter", arrayValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);

  RCP<const EnhancedNumberValidator< T > > readInPrototypeValidator =
    rcp_dynamic_cast<const ArrayValidator<EnhancedNumberValidator<T>, T > >(
      readInPL->getEntry(
        arrayParameterName).validator(), true)->getPrototype();
  RCP<const EnhancedNumberValidator< T > > actualPrototypeValidator =
    arrayValidator->getPrototype();

  TEST_EQUALITY(
    readInPrototypeValidator->getMin(),
    actualPrototypeValidator->getMin()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->getMax(),
    actualPrototypeValidator->getMax()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->getStep(),
    actualPrototypeValidator->getStep()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->getPrecision(),
    actualPrototypeValidator->getPrecision()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->hasMin(),
    actualPrototypeValidator->hasMin()
  );
  TEST_EQUALITY(
    readInPrototypeValidator->hasMax(),
    actualPrototypeValidator->hasMax()
  );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(Teuchos_Validator, StringToIntegralConverterTest, T)
{
  std::string defaultStringToIntegralParameterName = "defaultsti";
  std::string stringToIntegralParameterName = "sti";
  ParameterList myList;
  RCP<StringToIntegralParameterEntryValidator< T > > defaultStiValidator = rcp(
    new StringToIntegralParameterEntryValidator< T >(
      tuple<std::string>("value1", "value2", "value3"), stringToIntegralParameterName));
  RCP<StringToIntegralParameterEntryValidator< T > > stiValidator = rcp(
    new StringToIntegralParameterEntryValidator< T >(
      tuple<std::string>("value3", "value4", "value5"), 
      tuple<std::string>("the third value", "the fourth value", "the fifth value"),
      tuple< T >(3,4,5),
      stringToIntegralParameterName));
  myList.set(defaultStringToIntegralParameterName,
    "value1", "parameter with default sti validator", defaultStiValidator);
  myList.set(stringToIntegralParameterName, "value3", "parameter with sti validator",
    stiValidator);

  RCP<ParameterList> readInPL = writeThenReadPL(myList);


  RCP<const StringToIntegralParameterEntryValidator< T > > 
  readInDefaultStiValidator =
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator< T > >(
      readInPL->getEntry(
        defaultStringToIntegralParameterName).validator(), true);
  RCP<const StringToIntegralParameterEntryValidator< T > > 
  readInStiValidator =
    rcp_dynamic_cast<const StringToIntegralParameterEntryValidator< T > >(
      readInPL->getEntry(
        stringToIntegralParameterName).validator(), true);

  Array<std::string> readInDefaultValidStrings =
    *(readInDefaultStiValidator->validStringValues());
  Array<std::string> defaultValidStrings =
    *(defaultStiValidator->validStringValues());
  TEST_COMPARE_ARRAYS(readInDefaultValidStrings, defaultValidStrings);

  TEST_ASSERT(readInDefaultStiValidator->getStringDocs().is_null());
  TEST_EQUALITY( readInDefaultStiValidator->getDefaultParameterName(),
    defaultStiValidator->getDefaultParameterName());
  for(int i=0; i<defaultValidStrings.size(); ++i){
    TEST_EQUALITY(defaultStiValidator->getIntegralValue(defaultValidStrings[i]),
      readInDefaultStiValidator->getIntegralValue(defaultValidStrings[i]));
  }

  Array<std::string> readInValidStrings = *(readInStiValidator->validStringValues());
  Array<std::string> validStrings = *(stiValidator->validStringValues());
  TEST_COMPARE_ARRAYS(readInValidStrings, validStrings);

  TEST_COMPARE_ARRAYS(*(readInStiValidator->getStringDocs()),
    *(stiValidator->getStringDocs()));
  TEST_EQUALITY( readInStiValidator->getDefaultParameterName(),
    stiValidator->getDefaultParameterName());
  for(int i=0; i<validStrings.size(); ++i){
    TEST_EQUALITY(stiValidator->getIntegralValue(validStrings[i]),
      readInStiValidator->getIntegralValue(validStrings[i]));
  }

}

TEUCHOS_UNIT_TEST(Teuchos_Validator, existingPrototypeTest){
  ParameterList pl("ExsitingPrototypeList");
  RCP<StringValidator> stringVali = rcp(new StringValidator());
  RCP<ArrayValidator<StringValidator, std::string> > arrayStringVali 
    = rcp(new ArrayValidator<StringValidator, std::string>(stringVali));
  Array<std::string> strArray = tuple<std::string>("blah", "blah", "blah");
  pl.set("string param", "hi", "a string param", stringVali);
  pl.set("string array param", strArray, 
    "a string array parameter", arrayStringVali);
  RCP<ParameterList> readInPL = writeThenReadPL(pl);
  RCP<const ArrayValidator<StringValidator, std::string> > 
    inArrayValidator = 
    rcp_dynamic_cast<const ArrayValidator<StringValidator, std::string> >(
      readInPL->getEntry("string array param").validator(), true);
  TEST_ASSERT(readInPL->getEntry("string param").validator() 
    == inArrayValidator->getPrototype());
}


#define FULL_NUMBER_TYPE_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, EnhancedNumberValidatorConverter, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, NumberArrayValidatorConverterTest, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, StringToIntegralConverterTest, T ) 

#define NONINTEGRAL_NUMBER_TYPE_TEST( T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, EnhancedNumberValidatorConverter, T ) \
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(Teuchos_Validator, NumberArrayValidatorConverterTest, T ) 

typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned long ulong;


FULL_NUMBER_TYPE_TEST(int)
FULL_NUMBER_TYPE_TEST(uint)
FULL_NUMBER_TYPE_TEST(short)
FULL_NUMBER_TYPE_TEST(ushort)
FULL_NUMBER_TYPE_TEST(long)
FULL_NUMBER_TYPE_TEST(ulong)
NONINTEGRAL_NUMBER_TYPE_TEST(double)
NONINTEGRAL_NUMBER_TYPE_TEST(float)
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int llint;
typedef unsigned long long int ullint;
FULL_NUMBER_TYPE_TEST(llint)
FULL_NUMBER_TYPE_TEST(ullint)
#endif


} // namespace Teuchos

