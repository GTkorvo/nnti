//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $Zoltan_QueryContainer.h$
//
// Purpose        : Static Container object to allow Static (C-style)
//                  functions to access methods of dynamic objects.
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra
//
// Creation Date  : 08/04/2000
//
// Revision Information:
// ---------------------
//
// Revision Number: $$
//
// Revision Date  : $$
//
// Current Owner  : $$
//-------------------------------------------------------------------------

#ifndef ZOLTAN_QUERYCONTAINER_H_
#define ZOLTAN_QUERYCONTAINER_H_

#include <map>

class Zoltan_QueryObject;

class Zoltan_QueryContainer
{

public:

  static void setQueryID( const int & id );

  static const int & getQueryID();

  static bool registerQueryObject( const int & id,
				Zoltan_QueryObject * obj_ptr );

  static Zoltan_QueryObject * getQueryObject( const int & id );

private:

  static int CurrentObject;

  static std::map< int, Zoltan_QueryObject * > StaticMap;

};

#endif
