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

#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"


using namespace Teuchos;


Array<std::string> StrUtils::readFile(std::istream& is, char comment)
{
  std::string line;
  Array<std::string> rtn(0);

  while (readLine(is, line))
    {
      if (line.length() > 0) rtn.append(before(line, comment));
      line="";
    }
	
  return rtn;
}

Array<std::string> StrUtils::splitIntoLines(const std::string& input)
{
  int begin = 0;
  Array<std::string> rtn;
  const unsigned int len = input.length();
  for (unsigned int p=0; p<len; ++p) {
    const bool isEnd = p==len-1;
    if( input[p]=='\n' || input[p]=='\0' || input[p]=='\r' || isEnd )
    {
      if (p-begin > 1)
        rtn.append(
          subString( input, begin, p+(isEnd?(input[len-1]=='\n'?0:1):0) )
          );
      begin = p+1;
    }
  }
  return rtn;
}

Array<Array<std::string> > StrUtils::tokenizeFile(std::istream& is, char comment)
{
  std::string line;
  Array<Array<std::string> > rtn(0);
  Array<std::string> lines = readFile(is, comment);
  rtn.reserve(lines.length());
	
  int count = 0;
  for (int i=0; i<lines.length(); i++)
    {
      if (lines[i].length() == 0) continue;
      Array<std::string> tokens = stringTokenizer(lines[i]);
      if (tokens.length() == 0) continue;
      rtn.append(tokens);
      count++;
    }
	
  return rtn;
}

bool StrUtils::readLine(std::istream& is, std::string& line)
{
  char c[500];
  if (line.length() > 0) line[0] = '\0';
	
  if (is.eof()) return false;
  if (is.getline(c, 499))
    {
      line = std::string(c);
    }
	
  return true;
}

	

Array<std::string> StrUtils::getTokensPlusWhitespace(const std::string& str){
  Array<std::string> rtn(0);
  unsigned int start = 0;
	
  while(start < str.length())
    {
      unsigned int wordStart =  findNextNonWhitespace(str, start);
			/* add any preceding whitespace */
			if (wordStart > start)
				{
					rtn.append(subString(str, start, wordStart));
				}
			start = wordStart;
			/* add the next word */
      int stop = findNextWhitespace(str, start);
      if (start-stop == 0) return rtn;
      std::string sub = subString(str, start, stop);
      rtn.append(sub);
			start = stop;// findNextNonWhitespace(str, stop);
    }
  return rtn;
}

Array<std::string> StrUtils::stringTokenizer(const std::string& str){
  Array<std::string> rtn(0);
  unsigned int start = 0;
	
  while(start < str.length())
    {
      start =  findNextNonWhitespace(str, start);
      int stop = findNextWhitespace(str, start);
      if (start-stop == 0) return rtn;
      std::string sub = subString(str, start, stop);
      rtn.append(sub);
      start =  findNextNonWhitespace(str, stop);
    }
  return rtn;
}

std::string StrUtils::reassembleFromTokens(const Array<std::string>& tokens, int iStart)
{
  std::string rtn;

  for (int i=iStart; i<tokens.length(); i++) 
    {
      rtn += tokens[i];
      if (i < (tokens.length()-1)) rtn += " ";
    }
  return rtn;
}

void StrUtils::splitList(const std::string& big, Array<std::string>& list) 
{
  if (subString(big, 0,1)!="[") 
    {
      list.resize(1);
      list[0] = big;
      return;
    }
	
  int parenDepth = 0;
  int localCount = 0;
  std::string tmp(big);
  list.resize(0);

  // start at 1 to ignore '[';
	
  for (unsigned int i=1; i<big.length(); i++)
    {
      if (big[i]=='(') parenDepth++;
      if (big[i]==')') parenDepth--;
      if (big[i]==']') 
	{
	  tmp[localCount]='\0'; 
	  list.append(tmp);
	  break;
	}
      if (big[i]==',' && parenDepth==0)
	{
	  tmp[localCount]='\0';
	  list.append(tmp);
	  tmp = big;
	  localCount = 0;
	  continue;
	}
      tmp[localCount] = big[i];
      localCount++;
    }
}
							

// return the position of the next whitespace in a std::string. 
// If no whitespace, return -1;

int StrUtils::findNextWhitespace(const std::string& str, int offset)
{
  for (unsigned int i=0; i<(str.length()-offset); i++)
    {
      if (str[i+offset]==' ' || str[i+offset]=='\t' || str[i+offset]=='\n')
	{
	  return i+offset;
	}
    }
  return str.length();
}

int StrUtils::findNextNonWhitespace(const std::string& str, int offset)
{
  for (unsigned int i=0; i<(str.length()-offset); i++)
    {
      if (!(str[i+offset]==' ' || str[i+offset]=='\t' || str[i+offset]=='\n'))
	{
	  return i+offset;
	}
    }
  return str.length();
}


std::string StrUtils::varTableSubstitute(const std::string& rawLine,
				    const Array<std::string>& varNames,
				    const Array<std::string>& varValues)
{
  TEST_FOR_EXCEPTION(varNames.length() != varValues.length(),
                     std::runtime_error,
                     "mismatched variable tables in varTableSubstitute");
                     
  std::string line = rawLine;
  for (int i=0; i<varNames.length(); i++)
    {
      line = varSubstitute(line, varNames[i], varValues[i]);
    }
  return line;
}




std::string StrUtils::varSubstitute(const std::string& rawLine, 
			       const std::string& varName, 
			       const std::string& varValue)
{
  std::string line = rawLine;
  
  // iterate because there might be more than one occurance on this line
  while (find(line, varName) >= 0)
    {
      std::string b = before(line, varName);
      std::string a = after(line, varName);
      line = b + varValue + a;
    }
  return line;
}


std::string StrUtils::before(const std::string& str, char sub)
{
  char c[2];
  c[0] = sub;
  c[1] = 0;
  return before(str, c);
}

std::string StrUtils::before(const std::string& str, const std::string& sub)
{
  TEST_FOR_EXCEPTION(sub.c_str()==0,
                     std::runtime_error, "String::before: arg is null pointer");

  char* p = std::strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return str;
  int subLen = p-str.c_str();
  std::string rtn(str.c_str(), subLen);
  return rtn;
}

std::string StrUtils::after(const std::string& str, const std::string& sub)
{
  TEST_FOR_EXCEPTION(sub.c_str()==0,
                     std::runtime_error, "String::after: arg is null pointer");

  // find beginning of substring
  char* p = std::strstr((char*) str.c_str(), (char*) sub.c_str()) ;
  // if substring not found, return empty std::string
  if (p==0) return std::string();
  // offset to end of substring
  p+= std::strlen(sub.c_str());
  return std::string(p);
}

int StrUtils::find(const std::string& str, const std::string& sub)
{
  char* p = std::strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return -1;
  return p-str.c_str();
}

bool StrUtils::isWhite(const std::string& str)
{
  for (unsigned int i=0; i<str.length(); i++)
    {
      unsigned char c = str[i];
      if (c >= 33 && c <= 126)
        {
          return false;
        }
    }
  return true;
}

std::string StrUtils::fixUnprintableCharacters(const std::string& str)
{
  std::string rtn = str;
  for (unsigned int i=0; i<rtn.length(); i++)
    {
      unsigned char c = rtn[i];
      if (c < 33 || c > 126) 
        {
          if (c != '\t' && c != '\n'&& c != '\r' && c != '\f' && c != ' ')
            {
              rtn[i] = ' ';
            }
        }
    }
  return rtn;
}

std::string StrUtils::between(const std::string& str, const std::string& begin,
			 const std::string& end, std::string& front,
			 std::string& back)
{
  front = before(str, begin);
  std::string middle = before(after(str, begin), end);
  back = after(str, end);
  return middle;
}


std::string StrUtils::subString(const std::string& str, int begin, int end)
{
return std::string(str.c_str()+begin, end-begin);
}

std::string StrUtils::readFromStream(std::istream& is)
{
  TEST_FOR_EXCEPTION(true, std::logic_error, 
                     "StrUtils::readFromStream isn't implemented yet");

	return "";
}

std::string StrUtils::allCaps(const std::string& s)
{
  std::string rtn = s;
  for (unsigned int i=0; i<rtn.length(); i++)
    {
      rtn[i] = toupper(rtn[i]);
    }
  return rtn;
}

double StrUtils::atof(const std::string& s)
{
	return std::atof(s.c_str());
}

int StrUtils::atoi(const std::string& s)
{
	return std::atoi(s.c_str());
}

std::ostream& StrUtils::printLines(
  std::ostream             &os
  ,const std::string       &linePrefix
  ,const std::string       &lines
  )
{
  typedef Teuchos::Array<std::string> array_t;
  array_t linesArray = splitIntoLines(lines);
  for( int i = 0; i < static_cast<int>(linesArray.size()); ++i )
  {
    os << linePrefix << linesArray[i] << "\n";
  }
  return os;
}
  

	
