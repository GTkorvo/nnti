/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

/*!
 *  \file ml_MultiLevelPreconditioner_XML.cpp
 *
 *  \brief Converter from an XML file to internally stored Teuchos::ParameterList.
 *
 *  \author Marzio Sala, ETHZ/D-INFK.
 *
 *  \date Last update to Doxygen: 19-Mar-06.
 *
 */

#include "ml_common.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_epetra.h"
#include "Epetra_Comm.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListReader.hpp"

using namespace Teuchos;

// ============================================================================ 
static void AddParameter(Teuchos::ParameterList& List, 
                         const string& name, 
                         const Teuchos::ParameterEntry& entry)
{
  if (entry.isType<int>())
  {
    List.set(name, getValue<int>(entry));
  }
  else if (entry.isType<double>())
  {
    List.set(name, getValue<double>(entry));
  }
  else if (entry.isType<float>())
  {
    List.set(name, getValue<float>(entry));
  }
  else if (entry.isType<string>())
  {
    List.set(name, getValue<string>(entry));
  }
  else if (entry.isType<char>())
  {
    List.set(name, getValue<char>(entry));
  }
  else if (entry.isType<bool>())
  {
    List.set(name, getValue<bool>(entry));
  }
}

// ============================================================================ 
static void AddSubList(Teuchos::ParameterList& List, Teuchos::ParameterList& ListToAdd)
{
  if (List.name() == "ANONYMOUS") List.setName(ListToAdd.name());
  for (ParameterList::ConstIterator i = ListToAdd.begin(); i != ListToAdd.end(); ++i)
  {
    const ParameterEntry& val = ListToAdd.entry(i);
    const string& name = ListToAdd.name(i);
    if (val.isList()) AddSubList(List.sublist(name), ListToAdd.sublist(name));
    AddParameter(List, name, val);
  }
}

// ============================================================================ 
int ML_Epetra::ReadXML(const string &FileName, ParameterList &List,
            const Epetra_Comm &Comm)
{
  int i = 0, j;
  FILE* ML_capture_flag;
  ML_capture_flag = fopen((char*)FileName.c_str(),"r");
  if(ML_capture_flag) 
  {
    i++;
    fclose(ML_capture_flag);
  }

  Comm.SumAll(&i, &j, 1);

  int OutputLevel=0;
  if  (List.isParameter("ML output")) OutputLevel = List.get("ML output", 0);
  else if (List.isParameter("output")) OutputLevel = List.get("output", 0);

  if (j == 0) {
    if (OutputLevel && Comm.MyPID() == 0)
      cout << "***" << endl
           << "*** Unable to open XML input file '" << FileName << "'" << endl
           << "***" << endl;
    return(0);
  }
  
  Teuchos::FileInputSource fileSrc(FileName);
  Teuchos::XMLObject fileXML = fileSrc.getObject();

  Teuchos::ParameterList ListToAdd;

  if (fileXML.getTag() == "ParameterList" &&
      fileXML.getRequired("name") == "MultiLevelPreconditioner")
  {
    Teuchos::XMLParameterListReader ListReader;
    ListToAdd = ListReader.toParameterList(fileXML);
  }

  int ao=-1; //append or overwrite
  if (ListToAdd.get("ResetList", false))
  {
    ao = 0;
    Teuchos::ParameterList NewList;
    List = NewList;
  }
  else
    ao = 1;

  string xxx = ListToAdd.get("SetDefaults", "not-set");
  if (xxx != "not-set") {
    ML_Epetra::SetDefaults(xxx, ListToAdd,0,0,false);
  }

  AddSubList(List, ListToAdd);

  if  (List.isParameter("ML output")) OutputLevel = List.get("ML output", 0);
  else if (List.isParameter("output")) OutputLevel = List.get("output", 0);
  if ( (5 < OutputLevel) && (Comm.MyPID() == 0) ) {
    cout << "***" << endl;
    cout << "***" << " Reading XML file `" << FileName << "'..." << endl;
    if (ao == 0)
      cout << "***" << " Reset stored list" << endl;
    else if (ao == 1)
      cout << "***" << " Parameters are added to the stored list" << endl;
    if (xxx != "not-set")
      cout << "***" << " Setting default values to type `"
           << xxx << "'" << endl;
    cout << "***" << endl;
  }

  return(1);
}

#endif
