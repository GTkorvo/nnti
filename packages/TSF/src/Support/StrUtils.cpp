#include "StrUtils.h"


using namespace TSF;


TSFArray<string> StrUtils::readFile(istream& is, char comment)
{
  string line;
  TSFArray<string> rtn(0);

  while (readLine(is, line))
    {
      if (line.length() > 0) rtn.append(before(line, comment));
      line="";
    }
	
  return rtn;
}

TSFArray<string> StrUtils::splitIntoLines(const string& input)
{
  int begin = 0;
  TSFArray<string> rtn;

  for (unsigned int p=0; p<input.length(); p++)
    {
      if (input[p]=='\n' || input[p]=='\0' || input[p]=='\r')
	{
	  if (p-begin > 1) rtn.append(subString(input, begin, p));
	  begin = p+1;
	}
    }
  return rtn;
}
	

TSFArray<TSFArray<string> > StrUtils::tokenizeFile(istream& is, char comment)
{
  string line;
  TSFArray<TSFArray<string> > rtn(0);
  TSFArray<string> lines = readFile(is, comment);
  rtn.reserve(lines.length());
	
  int count = 0;
  for (int i=0; i<lines.length(); i++)
    {
      if (lines[i].length() == 0) continue;
      TSFArray<string> tokens = stringTokenizer(lines[i]);
      if (tokens.length() == 0) continue;
      rtn.append(tokens);
      count++;
    }
	
  return rtn;
}

bool StrUtils::readLine(istream& is, string& line)
{
  char c[500];
  if (line.length() > 0) line[0] = '\0';
	
  if (is.eof()) return false;
  if (is.getline(c, 499))
    {
      line = string(c);
    }
	
  return true;
}

	

TSFArray<string> StrUtils::getTokensPlusWhitespace(const string& str){
  TSFArray<string> rtn(0);
  unsigned int start = 0;
	
  while(start < str.length())
    {
      int wordStart =  findNextNonWhitespace(str, start);
			/* add any preceding whitespace */
			if (wordStart > start)
				{
					rtn.append(subString(str, start, wordStart));
				}
			start = wordStart;
			/* add the next word */
      int stop = findNextWhitespace(str, start);
      if (start-stop == 0) return rtn;
      string sub = subString(str, start, stop);
      rtn.append(sub);
			start = stop;// findNextNonWhitespace(str, stop);
    }
  return rtn;
}

TSFArray<string> StrUtils::stringTokenizer(const string& str){
  TSFArray<string> rtn(0);
  unsigned int start = 0;
	
  while(start < str.length())
    {
      start =  findNextNonWhitespace(str, start);
      int stop = findNextWhitespace(str, start);
      if (start-stop == 0) return rtn;
      string sub = subString(str, start, stop);
      rtn.append(sub);
      start =  findNextNonWhitespace(str, stop);
    }
  return rtn;
}

string StrUtils::reassembleFromTokens(const TSFArray<string>& tokens, int iStart)
{
  string rtn;

  for (int i=iStart; i<tokens.length(); i++) 
    {
      rtn += tokens[i];
      if (i < (tokens.length()-1)) rtn += " ";
    }
  return rtn;
}

void StrUtils::splitList(const string& big, TSFArray<string>& list) 
{
  if (subString(big, 0,1)!="[") 
    {
      list.resize(1);
      list[0] = big;
      return;
    }
	
  int parenDepth = 0;
  int localCount = 0;
  string tmp(big);
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
							

// return the position of the next whitespace in a string. 
// If no whitespace, return -1;

int StrUtils::findNextWhitespace(const string& str, int offset)
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

int StrUtils::findNextNonWhitespace(const string& str, int offset)
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


string StrUtils::varTableSubstitute(const string& rawLine,
				    const TSFArray<string>& varNames,
				    const TSFArray<string>& varValues)
{
  if (varNames.length() != varValues.length())
    TSFError::raise("mismatched variable tables in varTableSubstitute");

  string line = rawLine;
  for (int i=0; i<varNames.length(); i++)
    {
      line = varSubstitute(line, varNames[i], varValues[i]);
    }
  return line;
}




string StrUtils::varSubstitute(const string& rawLine, 
			       const string& varName, 
			       const string& varValue)
{
  string line = rawLine;
  
  // iterate because there might be more than one occurance on this line
  while (find(line, varName) >= 0)
    {
      string b = before(line, varName);
      string a = after(line, varName);
      line = b + varValue + a;
    }
  return line;
}


string StrUtils::before(const string& str, char sub)
{
  char c[2];
  c[0] = sub;
  c[1] = 0;
  return before(str, c);
}

string StrUtils::before(const string& str, const string& sub)
{
  if (sub.c_str()==0) TSFError::raise("String::before: arg is null pointer");
  char* p = strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return str;
  int subLen = p-str.c_str();
  string rtn(str.c_str(), subLen);
  return rtn;
}

string StrUtils::after(const string& str, const string& sub)
{
  if (sub.c_str()==0) TSFError::raise("String::after: arg is null pointer");
  // find beginning of substring
  char* p = strstr((char*) str.c_str(), (char*) sub.c_str()) ;
  // if substring not found, return empty string
  if (p==0) return string();
  // offset to end of substring
  p+= strlen(sub.c_str());
  return string(p);
}

int StrUtils::find(const string& str, const string& sub)
{
  char* p = strstr((char*) str.c_str(), (char*) sub.c_str());
  if (p==0) return -1;
  return p-str.c_str();
}

bool StrUtils::isWhite(const string& str)
{
// look for non-white characters
  int nWhite = (int) strspn(str.c_str(), " \n\t\r");
  //int n = (int) strcspn(str.c_str(), " \n\t\r");
  //if (n==0) return true;
  if (nWhite < str.length()) return false;
  return true;
}

string StrUtils::between(const string& str, const string& begin,
			 const string& end, string& front,
			 string& back)
{
  front = before(str, begin);
  string middle = before(after(str, begin), end);
  back = after(str, end);
  return middle;
}


string StrUtils::subString(const string& str, int begin, int end)
{
return string(str.c_str()+begin, end-begin);
}

string StrUtils::readFromStream(istream& is)
{
	TSFError::raise("StrUtils::readFromStream isn't implemented yet");
	return "";
}

string StrUtils::allCaps(const string& s)
{
	string rtn = s;
	for (unsigned int i=0; i<rtn.length(); i++)
		{
			rtn[i] = toupper(rtn[i]);
		}
	return rtn;
}

double StrUtils::atof(const string& s)
{
	return ::atof(s.c_str());
}

int StrUtils::atoi(const string& s)
{
	return ::atoi(s.c_str());
}



	
