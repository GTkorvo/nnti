/*! \file misc.c
    \brief Core routines for CMMLIB.
    \ingroup base
*/

/*!
   \defgroup base
*/

/*!
   \defgroup blas
*/

#include "cmm.h"
#ifndef DOXYGEN_SKIP
#include "stdlib.h"
#include "stdarg.h"
#endif

/*!
    \brief Returns size in bytes of given type.
    \param type  See definitions in cmm.h
*/
int   cmm_typesize(int type)
  {
   switch(type)
     {
      case CMM_BYTE:           return sizeof(unsigned char);
      case CMM_CHAR:           return sizeof(char);
      case CMM_UINT:           return sizeof(int);
      case CMM_INTEGER:        return sizeof(int);
      case CMM_SINGLE_REAL:    return sizeof(float);
      case CMM_REAL:           return sizeof(r8);
      case CMM_SINGLE_COMPLEX: return sizeof(float)*2;
      case CMM_COMPLEX:        return sizeof(c8);
      default:
        return 0;
     }
  }

/*!
    \brief Returns name of given type.
    \param type  See definitions in cmm.h

    Should be const char *

*/
char *cmm_typename(int type)
  {
   switch(type)
     {
      case CMM_BYTE:           return "Byte";
      case CMM_CHAR:           return "Char";
      case CMM_UINT:           return "UInt";
      case CMM_INTEGER:        return "Int";
      case CMM_SINGLE_REAL:    return "r4";
      case CMM_REAL:           return "r8";
      case CMM_SINGLE_COMPLEX: return "c4";
      case CMM_COMPLEX:        return "c8";
      default:
        return NULL;
     }
  }

/*!
    \brief Terminate program with status.
    \param status Exit status assigned to program.
*/
void cmm_exit(int status)
{
  exit(status);
}

/*!
    \brief Print message on stdio.
    \param fmt Format string.
    \param ... Print arguments.

   (see documentation for printf).
*/
void cmm_print(char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
}

/*!
    \brief Print message on stderr.
    \param fmt Format string.
    \param ... Print arguments.

   (see documentation for printf).
*/
void cmm_warn(char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
}

/*!
    \brief Print message on stderr and exit.
    \param fmt Format string.
    \param ... Print arguments.

   (see documentation for printf).
*/
void cmm_fatal(char *fmt, ...)
{
  char line[1024];

  va_list ap;

  va_start(ap, fmt);
  vsprintf(line, fmt, ap);
  va_end(ap);

  fprintf(stderr, "FATAL ERROR: %s\nending program.\n\n",line);
  cmm_exit(1);
}

/*!
    \brief Print message on stdio and wait for user response.
    \param fmt Format string.
    \param ... Print arguments.

   (see documentation for printf).
*/
void cmm_pause(char *fmt, ...)
{
  char line[80];
  va_list ap;

  va_start(ap, fmt);
  fprintf(stderr,"\n");
  vfprintf(stderr, fmt, ap);
  fprintf(stderr,"\n");
  va_end(ap);

  printf("pause>");
  fgets(line, 80, stdin);

}

