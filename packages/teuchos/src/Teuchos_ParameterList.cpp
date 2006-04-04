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

//#define TEUCHOS_PARAMETER_LIST_SHOW_TRACE

#include "Teuchos_ParameterList.hpp"

#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
#include "Teuchos_VerboseObject.hpp"
#endif

/* NOTE: ASCI Red (TFLOP) does not support the i-> function for iterators 
 * in the STL.  Therefore when compiling for the TFLOP we must redefine the 
 * iterator from i-> to (*i). This slows things down on other platforms 
 * so we switch between the two when necessary.
 */

namespace {

std::string filterValueToString(const Teuchos::ParameterEntry& entry )
{
  return ( entry.isList() ? std::string("...") : toString(entry.getAny()) );
}

struct ListPlusValidList {
  Teuchos::ParameterList   *list;
  Teuchos::ParameterList   *validList;
  ListPlusValidList(
    Teuchos::ParameterList   *_list
    ,Teuchos::ParameterList  *_validList
    )
    :list(_list),validList(_validList)
    {}
};

} // namespace 

namespace Teuchos {

ParameterList::ParameterList()
  :name_("ANONYMOUS")
{}

ParameterList::ParameterList(const std::string &name)
  :name_(name)
{}

ParameterList::ParameterList(const ParameterList& source) 
{
  name_ = source.name_;
  params_ = source.params_;
}

ParameterList& ParameterList::operator=(const ParameterList& source) 
{
  if (&source == this)
    return *this;
  name_ = source.name_;
  params_ = source.params_;
  return *this;
}

ParameterList& ParameterList::setParameters(const ParameterList& source) 
{
  for( ConstIterator i = source.begin(); i != source.end(); ++i ) {
    const std::string     &name_i  = this->name(i);
    const ParameterEntry  &entry_i = this->entry(i);
    if(entry_i.isList()) {
      this->sublist(name_i).setParameters(getValue<ParameterList>(entry_i));
    }
    else {
      this->setEntry(name_i,entry_i);
    }
  }
  this->updateSubListNames();
  return *this;
}

ParameterList::~ParameterList() 
{}

void ParameterList::unused(ostream& os) const
{
  for (ConstIterator i = params_.begin(); i != params_.end(); ++i) {
    if (!(entry(i).isUsed())) {
      os << "WARNING: Parameter \"" << name(i) << "\" " << entry(i)
         << " is unused" << endl;
    }
  }
}

std::string ParameterList::currentParametersString() const
{
  std::ostringstream oss;
  oss << "{";
  ParameterList::ConstIterator itr;
  int i;
  for( itr = this->begin(), i = 0; itr != this->end(); ++itr, ++i ) {
    const std::string     &entryName   = this->name(itr);
    const ParameterEntry  &entry       = this->entry(itr);
    if(i) oss << ",";
    oss
      << "\""<<entryName<<"\":"<<entry.getAny().type().name()
      <<"="<<filterValueToString(entry);
  }
  oss << "}";
  return oss.str();
}

bool ParameterList::isSublist(const string& name) const
{
  ConstIterator i = params_.find(name);

  if (i != params_.end())
    return (entry(i).isList());

  return false;
}

bool ParameterList::isParameter(const string& name) const
{
  return (params_.find(name) != params_.end());
}

ParameterList& ParameterList::sublist(const string& name)
{
  // Find name in list, if it exists.
  Iterator i = params_.find(name);

  // If it does exist and is a list, return the list value.
  // Otherwise, throw an error.
  if (i != params_.end()) {
    TEST_FOR_EXCEPTION(
      !entry(i).isList(), std::runtime_error,
      " Parameter " << name << " is not a list, it is of type \""
      <<entry(i).getAny().type().name()<<"\"!" );
    return getValue<ParameterList>(entry(i));
  }

  // The list does not exist so create a new empty list and return its reference
  const ParameterList newSubList(this->name()+std::string("->")+name);
  return any_cast<ParameterList>(
    params_.insert(
      Map::value_type(name,ParameterEntry(newSubList,false,true))
      ).first->second.getAny()
    );
  // Note that above I am very careful to construct the parameter list entry
  // directly in the insertion call to avoid the creation of a tempory list
  // object.  This looks ugly but it should be fast.
}

const ParameterList& ParameterList::sublist(const string& name) const
{
  // Find name in list, if it exists.
  ConstIterator i = params_.find(name);

  // If it does not exist, throw an error
  TEST_FOR_EXCEPTION( i == params_.end(), std::runtime_error,
                      " Parameter " << name << " is not a valid list!" );

  // If it does exist and is a list, return the list value.
  TEST_FOR_EXCEPTION( !entry(i).isList(), std::runtime_error,
                      " Parameter " << name << " is not a list!" );
  return getValue<ParameterList>(entry(i));
}
  
ostream& ParameterList::print(ostream& os, int indent, bool showTypes) const
{
  if (params_.begin() == params_.end()) {
    for (int j = 0; j < indent; j ++)
      os << ' ';
    os << "[empty list]" << endl;
  }
  else { 
    // Print parameters first
    for (ConstIterator i = params_.begin(); i != params_.end(); ++i) 
    {
      const ParameterEntry &entry_i = entry(i);
      if(entry_i.isList())
        continue;
      for (int j = 0; j < indent; j ++)
        os << ' ';
      os << name(i);
      if(showTypes)
        os << " : " << entry_i.getAny().type().name();
      os << " = " << entry_i << endl;
    }
    // Print sublists second
    for (ConstIterator i = params_.begin(); i != params_.end(); ++i) 
    {
      const ParameterEntry &entry_i = entry(i);
      if(!entry_i.isList())
        continue;
      for (int j = 0; j < indent; j ++)
        os << ' ';
      os << name(i) << " -> " << endl;
      getValue<ParameterList>(entry_i).print(os, indent + 2, showTypes);
    }
  }
  return os;
}

ParameterList::ConstIterator ParameterList::begin() const
{
  return params_.begin();
}

ParameterList::ConstIterator ParameterList::end() const
{
  return params_.end();
}



#if defined(TFLOP)

const string& ParameterList::name(ConstIterator i) const
{
  return ((*i).first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return ((*i).second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return ((*i).second);
}

#else

const string& ParameterList::name(ConstIterator i) const
{
  return (i->first);
}

ParameterEntry& ParameterList::entry(Iterator i)
{
  return (i->second);
}

const ParameterEntry& ParameterList::entry(ConstIterator i) const
{
  return (i->second);
}

#endif

void ParameterList::validateParameters(
  const ParameterList              &validParamList
  ,const int                       depth
  ,const EValidateUsed             validateUsed
  ,const EValidateDefaults         validateDefaults
  ) const
{
  typedef std::deque<ListPlusValidList> sublist_list_t;
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
  RefCountPtr<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  OSTab tab(out);
  *out << "\n*** Entering ParameterList::validateParameters(...) for this->name()=\""<<this->name()<<"\"...\n";
#endif
  //
  // First loop through and validate the parameters at this level.
  //
  // Here we generate a list of sublists that we will search next
  //
  sublist_list_t sublist_list;
  ConstIterator itr;
  for( itr = this->begin(); itr != this->end(); ++itr ) {
    const std::string    &entryName   = this->name(itr);
    const ParameterEntry &entry       = this->entry(itr);
    if(
      ( entry.isUsed() && validateUsed!=VALIDATE_USED_ENABLED )
      ||
      ( entry.isDefault() && validateDefaults!=VALIDATE_DEFAULTS_ENABLED )
      )
    {
      continue;
    }
    const ParameterEntry *validEntry = validParamList.getEntryPtr(entryName);
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
    OSTab tab(out);
    *out << "\nentryName=\""<<entryName<<"\"\n";
#endif
    const bool validType = ( validEntry!=NULL ? entry.getAny().type() == validEntry->getAny().type() : false );
    TEST_FOR_EXCEPTION(
      !( validEntry!=NULL && validType )
      ,Exceptions::InvalidParameter
      ,"Error, the parameter {name=\""<<entryName<<"\",type=\""<<entry.getAny().type().name()<<"\""
      ",value=\""<<filterValueToString(entry)<<"\"} in the parameter (sub)list \""<<this->name()<<"\" "
      << ( validEntry==NULL
           ? "was not found in the list of valid parameters"
           : std::string("exists in the list of valid parameters but has the wrong type.  The correct type is \"")
           +validEntry->getAny().type().name()+std::string("\"")
        )
      << ". The valid parameters and types are "<<validParamList.currentParametersString()
      );
    if( entry.isList() && depth > 0 ) {
      sublist_list.push_back(
        ListPlusValidList(
          &getValue<ParameterList>(entry),&getValue<ParameterList>(*validEntry)
          )
        );
    }
  }
  //
  // Now loop through the sublists and validate their parameters
  //
  for( sublist_list_t::const_iterator sl_itr = sublist_list.begin(); sl_itr != sublist_list.end(); ++sl_itr ) {
    sl_itr->list->validateParameters(
      *sl_itr->validList
      ,depth-1
      ,validateUsed
      ,validateDefaults
      );
  }
#ifdef TEUCHOS_PARAMETER_LIST_SHOW_TRACE
  *out << "\n*** Existing ParameterList::validateParameters(...) for this->name()=\""<<this->name()<<"\"...\n";
#endif
}

// private

void ParameterList::updateSubListNames(int depth)
{
  const std::string this_name = this->name();
  Map::iterator itr;
  for( itr = params_.begin(); itr != params_.end(); ++itr ) {
    const std::string    &entryName   = this->name(itr);
    const ParameterEntry &entry       = this->entry(itr);
    if(entry.isList()) {
      ParameterList &sublist = getValue<ParameterList>(entry);
      sublist.setName(this_name+std::string("->")+entryName);
      if(depth > 0)
        sublist.updateSubListNames(depth-1);
    }
  }
}

} // namespace Teuchos
