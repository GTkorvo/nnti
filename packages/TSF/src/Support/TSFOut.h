#ifndef TSFOUT_H
#define TSFOUT_H

#include "TSFDefs.h"
#include <string>
#include <cstdarg>
#include <stdexcept>
#include "TSFWriterBase.h"
#include "TSFSmartPtr.h"

namespace TSF
{

  using std::string;

  /** \ingroup IO
   *
   */
  class TSFOut
    {
    public:
      /** */
      static void print(const string& msg);

      /** */
      static void println(const string& msg);

      /** */
      static void rootPrintln(const string& msg);

      /** */
      static void printf(const char* format, ...);

      /** */
      static void rootPrintf(const char* format, ...);

      /** */
      static void setWriter(const TSFSmartPtr<TSFWriterBase>& writer );

    private:
      static TSFSmartPtr<TSFWriterBase> writer_;

      /** */
      static void vprintf(const char* format, va_list args);

			/** */
			static int hack_vsnprintf(char* str, size_t size, const char* format, va_list args);
    };
}

#endif
