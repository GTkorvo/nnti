/*! \file str.c
    \brief functions for string manipulation.
    \ingroup base
*/
#include "cmm.h"

#ifndef DOXYGEN_SKIP
#include <ctype.h>
#include <string.h>
#endif

/*!
   \brief Clone a string.
   \param s Source string.
   \return Cloned string.

    Allocates storage and then copies the string
    into it, returning the pointer to the new string.

    The allocated string should be free'd with cmm_free().

*/
char *cmm_str_clone(char *s)
{
  char *t;
  int l;

  if (s == 0)
    return cmm_str_clone("");

  l = strlen(s);
  t = (char *) cmm_alloc(l + 1);
  strcpy(t, s);
  return t;
}

/*!
  \brief Copy a string.
  \param src  Source string.
  \param dst  Destination string.

   \return The count of characters copied.

  Copies string from source buffer
  to destination buffer.

*/
int cmm_str_copy(char *src, char *dst)
{
  int l;

  if(src == NULL)
    {
     *dst = 0;
     return 0;
    }

  l = strlen(src);
  strcpy(dst, src);
  return l;
}

/*!
  \brief Copy a string up to n characters.
  \param src Source string.
  \param dst Destination string.
  \param n   Maximum number of characters to copy.
   \return The count of characters copied.

    Copies a string from source buffer \texttt{src} to
    to destination buffer \texttt{dst}. At most \texttt{n}
    characters are copied.


*/

int cmm_str_copyn(char *src, char *dst, int n)
    {
     int i, count;
     char ch;

  if(src == NULL)
    {
     *dst = 0;
     return 0;
    }

     count = 0;
     for(i=0;i<n;i++)
       {
        ch = *src++;
        if(ch == 0)
          break;
        *dst++ = ch;
        count++;
       }
     *dst = 0;
     return count+1;
    }

/*!
   \brief Copy a string up to delimiter.
   \param src Source string.
   \param dst Destination string.
   \param delims String of delimiter characters.
   \return Number of characters copied.

    Copy a string from source buffer \texttt{src} to
    to destination buffer \texttt{dst}. When any character
    in the delimiter string \texttt{delims} is encountered
    the copy is terminated.

*/
int cmm_str_copyd(char *src, char *dst, char *delims)
{
  int nd = strlen(delims);
  int l = strlen(src);
  int i, k;

  k = 0;
  for (i = 0; i < l; i++)
    {
      int j;
      int ch = *src++;

      if (ch == 0)
        break;
      for (j = 0; j < nd; j++)
        if (ch == delims[j])
          {
            ch = 0;
            break;
          }
      if (ch == 0)
        break;
      *dst++ = ch;
      k++;
    }
  *dst = 0;

  return k;

}

/*!
  \brief Compute string length.
  \param src String.
  \return Number of characters in string.

  Return value does not include terminating 0.

*/
int cmm_str_len(char *src)
{
  return strlen(src);
}

/*!
  \brief Trim end of line from string.
  \param src String.
  \return Number of characters preceding end of line.

*/
int cmm_str_trim_eol(char *src)
{
  int l = strlen(src);

/* remove possible \r\n */
  if (l > 0)
    {
      char ch;

      ch = src[l - 1];
      if (ch == '\r' || ch == '\n')
        l--;
    }

/* remove possible \r\n */
  if (l > 0)
    {
      char ch;

      ch = src[l - 1];
      if (ch == '\r' || ch == '\n')
        l--;
    }
  src[l] = 0;
  return l;
}

/*!
  \brief Compare two strings.
  \param s String 1
  \param t String 2

  \return True if the strings are identical, False otherwise.

*/
Boolean cmm_str_compare(char *s, char *t)
{

  while (1)
    {
      char chs, cht;

      chs = *s++;
      cht = *t++;
      if (chs != cht)
        {
          return False;
        }

      if (chs == 0)
        {
          return True;
        }
    }

}

/*!
  \brief Compare two strings.
  \param s String 1.
  \param t String 2.
  \param maxlen Maximum number of characters to compare.
  \return True if strings are identical up to \c maxlen characters.
*/
Boolean cmm_str_compareto(char *s, char *t, int maxlen)
{
  int i;
  char ch;

  i = 0;
  while (1)
    {
      ch = *s++;
      if (ch != *t++)
        return False;
      if (ch == 0)
        break;
      i++;
      if (i >= maxlen)
        break;
    }

  return True;
}

/*!
  \brief Convert string to real.
  \param s String to convert.
  \return Converted value.

*/
r8 cmm_str_cvt_r8(char *s)
{
  r8 v;

  v = 0.;
  if (s)
    sscanf(s, "%le", &v);
  return v;
}

/*!
   \brief Convert string to integer.
   \param s String to convert.
   \return Converted value.

*/
int cmm_str_cvt_int(char *s)
    {
     int iv = 0;
     sscanf(s, "%d", &iv);
     return iv;
    }

/*!
   \brief Convert string to integer given choices.
   \param str String to convert.
   \param nchoices Number of choice values.
   \param choices List of choice values.
   \return Choice selected.

   This function determines the match for the string
   among a list of choices.

*/
int cmm_str_cvt_choice(char *str, int nchoices, char **choices)
    {
     int i;

     for(i=0;i<nchoices;i++)
       {
        if(cmm_str_compare(str, choices[i]))
          return i;
       }

     return -1;  // none of the above
    }


/*!
   \brief Convert string to upper case.
   \param s String to convert.

   String conversion is in place.
*/
void cmm_str_toupper(char *s)
{
  int i;
  char ch;
  int l = strlen(s);

  for (i = 0; i < l; i++)
    {
      ch = s[i];
      ch = toupper(ch);
      s[i] = ch;
    }
}

/*!
   \brief Convert string to lower case.
   \param s String to convert.

   String conversion is in place.
*/
void cmm_str_tolower(char *s)
{
  int i;
  char ch;
  int l = strlen(s);

  for (i = 0; i < l; i++)
    {
      ch = s[i];
      ch = tolower(ch);
      s[i] = ch;
    }
}

/* return 0 or length of next arg */
// This decodes quoted strings.
static int get_one_arg(char *line, char *argbuf, int *la)
{
  int l, larg;
  char ch;

  l = 0;
  larg = 0;

  while (1)
    {
      ch = *line++;
      l++;
      if (ch == ' ' || ch == '\n' || ch == '\t')
        continue;
      break;
    }

  if (ch == 0)
    goto done;

/* quoted argument */
  if (ch == '"')
    {
      while (1)
        {
          ch = *line++;
          l++;
          if (ch == 0)
            {
              *argbuf = 0;
              cmm_warn("quoted string not ended by \".\n");
              return l;
            }
          if (ch == '"')
            {
              *argbuf = 0;
              *la = larg;
              return l;
            }

          if (ch == '\\')  // handle escape sequences
            {
              ch = *line++;
              l++;
              if (ch == 'n')
                ch = '\n';
              if (ch == 't')
                ch = '\t';
            }

          *argbuf++ = ch;
          larg++;
        }
    }

/* non-quoted string */
  *argbuf++ = ch;
  larg++;
  while (1)
    {
      ch = *line++;
      l++;
      if (ch == 0 || ch == ' ' || ch == '\n' || ch == '\t')
        {
          *argbuf = 0;
          goto done;
        }
      *argbuf++ = ch;
      larg++;
    }

done:
  *la = larg;
  return l;
}

/*!
  \brief Parse a string into an argument list.
  \param str String to parse.
  \param maxargs Maximum number of arguments allowed.
  \param maxbuf  Maximum length of buffer.
  \param args    Returned list of arguments.
  \param buf     Buffer.
  \return Number of parsed arguments.

   This function parses a string into a list of whitespace separated argument
   strings. Arguments may be quoted, in which case the quoted string is decoded.

*/
int cmm_str_arglist(char *str, int maxargs, int maxbuf, char **args, char *buf)
{
  char argbuf[256];
  char *bp;
  int nargs;

  nargs = 0;
  bp = buf;
  while (1)
    {
      int l, la;

      la = 0;
      l = get_one_arg(str, argbuf, &la);
//      printf("arg %d larg %d\n",nargs,la);
      if (la == 0)
        break;
      str += l;
      args[nargs++] = bp;
      cmm_str_copy(argbuf, bp);
      bp += la + 1;
    }

  args[nargs] = 0;
  return nargs;
}

/*!
   \brief Parse complete file name into path, file, extension.

   \param name File name to parse.
   \param path Path name.
   \param file File name without path or extension.
   \param ext  Extension.

   \return parse status.

   This routine parses a complete file name or path name. The
   \c path field is returned empty if there is no path part. Note that
   trailing / will be stripped, also // is treated as /.

   The \c file field will be returned empty if there is no file, and
   similarly for the \c ext field. Note that the extension includes
   the "." character.

*/
int cmm_str_parse_file(char *name, char *path, char *file, char *ext)
{
  int i, l, k1, k2;

  l = strlen(name);

/* Find location of last / in string */
  k1 = -1;
  for (i = 0; i < l; i++)
    {
      if (name[i] == '/')
        k1 = i;
    }

  if (k1 == -1)
    {
      *path = 0;
    }
  else
    {
      cmm_str_copyn(name, path, k1);
    }

/* Find first location of . after path */
  k2 = k1 + 1;
  for (i = k1 + 1; i < l; i++)
    {
      if (name[i] == '.')
        {
          k2 = i;
          break;
        }
      if (i + 1 == l)
        k2 = l;
    }
  cmm_str_copyn(name + k1 + 1, file, k2 - k1 - 1);

/* Get extension */
  cmm_str_copyn(name + k2, ext, l - k2);
  return 0;
}

/*!
   \brief  Encode a string into a quoted string.
   \param src String to encode.
   \param dst Output string.
   \param maxlen Size of output string.

   The output string begins and ends with quotes.

   Inside the string the following mappings take place:
LIST
    "   -> \"
    tab -> \t
    cr  -> \r
    eol -> \n
LIST.
*/
void cmm_str_encode_qs(char *src, char *dst, int maxlen)
    {
     char ch;

     maxlen -= 3; // for opening and closing quotes and the null.

     *dst++ = '"';
     while(maxlen > 0)
       {
        ch = *src++;

        if(ch == 0)
          break;

        if(ch == '"')
          {
           *dst++ = '\\';
           *dst++ = '"';
           maxlen--;
          }
        else if(ch == '\t')
          {
           *dst++ = '\\';
           *dst++ = 't';
           maxlen--;
          }
        else if(ch == '\r')
          {
           *dst++ = '\\';
           *dst++ = 'r';
           maxlen--;
          }
        else if(ch == '\n')
          {
           *dst++ = '\\';
           *dst++ = 'n';
           maxlen--;
          }
        else
          *dst++ = ch;
        maxlen--;
       }
     
     *dst++ = '"';
     *dst++ = 0;

    }

/*!
  \brief Copy memory.
  \param l Number of bytes to copy.
  \param src Source address.
  \param dst Destination address.

*/
void cmm_mem_copy(int l, void *src, void *dst)
    {
     memcpy(dst, src, l);
    }

