#ifndef EXPATHANDLERADAPTER_H
#define EXPATHANDLERADAPTER_H

#include "TSFDefs.h"

#ifdef HAVE_EXPAT
#include "TreeBuildingXMLHandler.h"
#include "TSFSmartPtr.h"

#include "expat.h"

extern "C"
{
  /** \ingroup XML Expat callback for start of an XML element. */
  void expatStartElementHandler(void* handler,
                                const XML_Char* name,
                                const XML_Char** attr);

  /** \ingroup XML Expat callback for end of an XML element. */
  void expatEndElementHandler(void* handler,
                              const XML_Char* name);

  /** \ingroup XML Expat callback for character data. */
  void expatCharacterDataHandler(void* handler,
                                 const XML_Char* s,
                                 int len);
};

#endif


#endif
