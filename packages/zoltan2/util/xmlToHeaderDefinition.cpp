// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file xmlToHeaderDefinition.cpp
    \brief Builds the parameter header file required at compile time.
 */

#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <cstring>
#include <sstream>

using std::cout;
using std::endl;
using std::string;
using std::ofstream;
using std::ifstream;
using std::ostringstream;

static string fixLine(char *inLine)
{
  string oldLine(inLine);
  string quote("\"");
  
  size_t pos = oldLine.find(quote);
  if (pos == string::npos)
    return oldLine;

  string newLine;
  string quotedQuote("\\\"");

  for (int i=0; i < oldLine.size(); i++){
    char c = oldLine[i];
    if (c == quote[0])
      newLine.append(quotedQuote);
    else
      newLine.append(&c,1);
  }
  return newLine;
}

int main(int argc, char  *argv[])
{
  if (argc < 3){
    cout << "Usage: " << argv[0] << " xmlFileName hppFileName" << endl;
    return 1;
  }

  char *xmlFile = argv[1];
  char *hppFile = argv[2];

  ofstream oFile;
  ifstream iFile;

  try{
    oFile.open(hppFile);
  }
  catch(std::exception &e){
    cout << "Error: " << e.what() << " " << hppFile << endl;
    return 1;
  }

  try{
    iFile.open(xmlFile);
  }
  catch(std::exception &e){
    oFile.close();
    cout << "Error: " << e.what() << " " << xmlFile << endl;
    return 1;
  }

  oFile << "// \n";
  oFile << "// This file was automatically generated by CMake\n";
  oFile << "// with the following command:\n";
  oFile << "// " << argv[0] << " " << argv[1] << " " << argv[2] << ".\n";
  oFile << "// \n";

  oFile << "#ifndef ZOLTAN2_PARAMETER_DEFINITION_HEADER\n";
  oFile << "#define ZOLTAN2_PARAMETER_DEFINITION_HEADER\n";
  oFile << "\n#define ZOLTAN2_XML_PARAMETER_STRING \"";

  char lineBuf[1024];
  iFile.clear();
  bool go=false;
  int sanity=10000;

  while (sanity--){
    iFile.getline(lineBuf, 1024);
    if (!iFile.good())
      break;
    if (go && (strlen(lineBuf) > 1)) {
      oFile << " \\\n  " << fixLine(lineBuf);
    }
    else if (!go){
      string line(lineBuf);
      size_t pos = line.find("ParameterList");
      if (pos == string::npos)
        continue;
      oFile << " \\\n  " << fixLine(lineBuf);
      go = true;
    }
  }

  if (go == false){
    iFile.close();
    oFile.close();
    cout << "Error: ParameterList XML definition not found." << endl;
    return 1;
  }
      
  iFile.close();

  oFile << "\"\n\n#endif  //ZOLTAN2_PARAMETER_DEFINITION_HEADER\n";

  oFile.close();

  return 0;
}

