// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include "MueLu_Exceptions.hpp"
#include "MueLu_MasterList.hpp"

namespace MueLu {

  Teuchos::RCP<const Teuchos::ParameterList> MasterList::List() {
    if (masterList_.is_null()) {
      masterList_ = Teuchos::getParametersFromXmlString(stringList_);
    }

    return masterList_;
  }

  Teuchos::RCP<Teuchos::ParameterList> MasterList::GetProblemSpecificList(std::string const & problemType) {

    if ( (problemType != problemType_) || problemSpecificList_.is_null() ) {
      if (DefaultProblemTypeLists_.find(problemType) != DefaultProblemTypeLists_.end()) {
        problemType_ = problemType;
        problemSpecificList_ = Teuchos::getParametersFromXmlString(DefaultProblemTypeLists_[problemType]);
      } else {
        //TODO provide valid problem types
        TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Invalid problem type " << problemType << ".");
      }
    }
    return problemSpecificList_;
  }
 
   std::string MasterList::interpretParameterName(const std::string& name, const std::string& value) {

    // used to concatenate the return string
    std::stringstream ss;

    // put in short cuts here!

    if (name == "verbosity") { 
      std::string verb = "none";
      if (value == "\"0\"") verb = "none";
      if (value == "\"1\"" || value == "\"2\"" || value == "\"3\"") verb = "low";
      if (value == "\"4\"" || value == "\"5\"" || value == "\"6\"") verb = "medium";
      if (value == "\"7\"" || value == "\"8\"") verb = "high";
      if (value == "\"9\"") verb = "extreme";
      if (value == "\"10\"") verb = "test";  
      verb = "\"" + verb + "\"";
      ss << "<Parameter name=\"verbosity\" type=\"string\" value=" << verb << "/>"; 
      return ss.str(); 
    }
    
    if (name == "cycle type") {
      std::stringstream temp1; temp1 << "\"" << "MGV" << "\"";
      std::stringstream temp2; temp2 << "\"" << "MGV" << "\"";
      if (value == temp1.str() ) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"; return ss.str(); }
      else if (value == temp2.str()) { ss << "<Parameter name=\"cycle type\" type=\"string\" value=\"W\"/>"; return ss.str(); }
      else TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MasterList::interpretParameterName, Line " << __LINE__ << ". "
                                           << "The parameter " << value << " is not supported by MueLu.");
      return ss.str();
    }    

    // energy minimization is enabled
    if (name == "multigrid algorithm") {
      std::stringstream temp; temp << "\"" << "1" << "\"";
      if (value == temp.str() ) { ss << "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"pg\"/>"; return ss.str(); }
    }

    if (name == "repartition: enable") {
      std::stringstream temp1; temp1 << "\"" << "1" << "\"";
      if (value == temp1.str()) {
        RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        *out << "WARNING: repartitioning in MueLu is different to MLs. Please refer to the MueLu users Manual for more information." << std::endl;
      }
    }
    
    // put in auto-generated code here


    if (name == "number of equations") { ss << "<Parameter name=\"number of equations\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "max levels") { ss << "<Parameter name=\"max levels\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "problem: symmetric") { ss << "<Parameter name=\"problem: symmetric\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "aggregation: drop tol") { ss << "<Parameter name=\"aggregation: drop tol\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    if (name == "print initial parameters") { ss << "<Parameter name=\"print initial parameters\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "print unused parameters") { ss << "<Parameter name=\"print unused parameters\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "sa: damping factor") { ss << "<Parameter name=\"sa: damping factor\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    if (name == "sa: eigenvalue estimate num iterations") { ss << "<Parameter name=\"sa: eigenvalue estimate num iterations\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: enable") { ss << "<Parameter name=\"repartition: enable\" type=\"bool\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: start level") { ss << "<Parameter name=\"repartition: start level\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: min rows per proc") { ss << "<Parameter name=\"repartition: min rows per proc\" type=\"int\" value=" << value << "/>"; return ss.str(); }      
    if (name == "repartition: max imbalance") { ss << "<Parameter name=\"repartition: max imbalance\" type=\"double\" value=" << value << "/>"; return ss.str(); }      
    return "";
  }
 
  Teuchos::RCP<Teuchos::ParameterList> MasterList::masterList_ = Teuchos::null;
  Teuchos::RCP<Teuchos::ParameterList> MasterList::problemSpecificList_ = Teuchos::null;
  std::string                          MasterList::problemType_ = "unknown";
  const std::string                    MasterList::stringList_ =
"<ParameterList name=\"MueLu\">"
  "<Parameter name=\"problem: type\" type=\"string\" value=\"unknown\"/>"
  "<Parameter name=\"verbosity\" type=\"string\" value=\"high\"/>"
  "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"max levels\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"cycle type\" type=\"string\" value=\"V\"/>"
  "<Parameter name=\"problem: symmetric\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"parameterlist: syntax\" type=\"string\" value=\"muelu\"/>"
  "<Parameter name=\"smoother: pre or post\" type=\"string\" value=\"both\"/>"
  "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
  "<Parameter name=\"smoother: pre type\" type=\"string\" value=\"RELAXATION\"/>"
  "<Parameter name=\"smoother: post type\" type=\"string\" value=\"RELAXATION\"/>"
  "<ParameterList name=\"smoother: params\"/>"
  "<ParameterList name=\"smoother: pre params\"/>"
  "<ParameterList name=\"smoother: post params\"/>"
  "<Parameter name=\"smoother: overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"smoother: pre overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"smoother: post overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"coarse: max size\" type=\"int\" value=\"2000\"/>"
  "<Parameter name=\"coarse: type\" type=\"string\" value=\"SuperLU\"/>"
  "<ParameterList name=\"coarse: params\"/>"
  "<Parameter name=\"coarse: overlap\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: type\" type=\"string\" value=\"uncoupled\"/>"
  "<Parameter name=\"aggregation: ordering\" type=\"string\" value=\"natural\"/>"
  "<Parameter name=\"aggregation: drop scheme\" type=\"string\" value=\"classical\"/>"
  "<Parameter name=\"aggregation: drop tol\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"aggregation: min agg size\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"aggregation: max agg size\" type=\"int\" value=\"-1\"/>"
  "<Parameter name=\"aggregation: max selected neighbors\" type=\"int\" value=\"0\"/>"
  "<Parameter name=\"aggregation: Dirichlet threshold\" type=\"double\" value=\"0.0\"/>"
  "<Parameter name=\"aggregation: enable phase 1\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 2a\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 2b\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: enable phase 3\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"aggregation: preserve Dirichlet points\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: export visualization data\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"old\"/>"
  "<ParameterList name=\"export data\"/>"
  "<Parameter name=\"print initial parameters\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"print unused parameters\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"transpose: use implicit\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
  "<Parameter name=\"semicoarsen: coarsen rate\" type=\"int\" value=\"3\"/>"
  "<Parameter name=\"sa: damping factor\" type=\"double\" value=\"1.33\"/>"
  "<Parameter name=\"sa: use filtered matrix\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"sa: calculate eigenvalue estimate\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"sa: eigenvalue estimate num iterations\" type=\"int\" value=\"10\"/>"
  "<Parameter name=\"filtered matrix: use lumping\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"filtered matrix: reuse eigenvalue\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"filtered matrix: reuse graph\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"emin: iterative method\" type=\"string\" value=\"cg\"/>"
  "<Parameter name=\"emin: num iterations\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"emin: num reuse iterations\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"emin: pattern\" type=\"string\" value=\"AkPtent\"/>"
  "<Parameter name=\"emin: pattern order\" type=\"int\" value=\"1\"/>"
  "<Parameter name=\"repartition: enable\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: partitioner\" type=\"string\" value=\"zoltan2\"/>"
  "<ParameterList name=\"repartition: params\"/>"
  "<Parameter name=\"repartition: start level\" type=\"int\" value=\"2\"/>"
  "<Parameter name=\"repartition: min rows per proc\" type=\"int\" value=\"800\"/>"
  "<Parameter name=\"repartition: max imbalance\" type=\"double\" value=\"1.2\"/>"
  "<Parameter name=\"repartition: remap parts\" type=\"bool\" value=\"true\"/>"
  "<Parameter name=\"repartition: remap num values\" type=\"int\" value=\"4\"/>"
  "<Parameter name=\"repartition: print partition distribution\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"repartition: rebalance P and R\" type=\"bool\" value=\"false\"/>"
  "<Parameter name=\"reuse: type\" type=\"string\" value=\"none\"/>"
"</ParameterList>"
;
  std::map<std::string,std::string> MasterList::DefaultProblemTypeLists_ = DefaultProblemStrings<std::string,std::string>
("Poisson-2D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Poisson-3D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"1\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-2D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"2\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("Elasticity-3D",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"number of equations\" type=\"int\" value=\"3\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"CHEBYSHEV\"/>"
          
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"sa\"/>"
          
    "</ParameterList>"
  )
("MHD",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"SCHWARZ\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"schwarz: overlap level\" type=\"int\" value=\"1\"/>"
        
        "<Parameter name=\"schwarz: combine mode\" type=\"string\" value=\"Zero\"/>"
        
        "<Parameter name=\"schwarz: use reordering\" type=\"bool\" value=\"false\"/>"
        
        "<Parameter name=\"subdomain solver name\" type=\"string\" value=\"RILUK\"/>"
        
    "<ParameterList name=\"subdomain solver parameters\">"
    
        "<Parameter name=\"fact: iluk level-of-fill\" type=\"int\" value=\"0\"/>"
        
        "<Parameter name=\"fact: absolute threshold\" type=\"double\" value=\"0.\"/>"
        
        "<Parameter name=\"fact: relative threshold\" type=\"double\" value=\"1.\"/>"
        
        "<Parameter name=\"fact: relax value\" type=\"double\" value=\"0.\"/>"
        
    "</ParameterList>"
  
    "</ParameterList>"
  
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"transpose: use implicit\" type=\"bool\" value=\"true\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"unsmoothed\"/>"
          
    "</ParameterList>"
  )
("ConvectionDiffusion",

    "<ParameterList name=\"MueLu\">"
    
            "<Parameter name=\"problem: symmetric\" type=\"bool\" value=\"false\"/>"
          
            "<Parameter name=\"smoother: type\" type=\"string\" value=\"RELAXATION\"/>"
          
    "<ParameterList name=\"smoother: params\">"
    
        "<Parameter name=\"relaxation: type\" type=\"string\" value=\"Gauss-Seidel\"/>"
        
    "</ParameterList>"
  
            "<Parameter name=\"aggregation: mode\" type=\"string\" value=\"new\"/>"
          
            "<Parameter name=\"multigrid algorithm\" type=\"string\" value=\"pg\"/>"
          
            "<Parameter name=\"sa: use filtered matrix\" type=\"bool\" value=\"true\"/>"
          
    "</ParameterList>"
  )
;
  std::map<std::string,std::string> MasterList::ML2MueLuLists_ = DefaultProblemStrings<std::string,std::string>

         ("default values","problem: type")
      
         ("ML output","verbosity")
      
         ("PDE equations","number of equations")
      
         ("max levels","max levels")
      
         ("prec type","cycle type")
      
         ("problem: symmetric","problem: symmetric")
      
         ("parameterlist: syntax","parameterlist: syntax")
      
         ("smoother: pre or post","smoother: pre or post")
      
         ("smoother: type","smoother: type")
      
         ("smoother: pre type","smoother: pre type")
      
         ("smoother: post type","smoother: post type")
      
         ("smoother: params","smoother: params")
      
         ("smoother: pre params","smoother: pre params")
      
         ("smoother: post params","smoother: post params")
      
         ("smoother: overlap","smoother: overlap")
      
         ("smoother: pre overlap","smoother: pre overlap")
      
         ("smoother: post overlap","smoother: post overlap")
      
         ("max size","coarse: max size")
      
         ("coarse: type","coarse: type")
      
         ("coarse: params","coarse: params")
      
         ("coarse: overlap","coarse: overlap")
      
         ("aggregation: type","aggregation: type")
      
         ("aggregation: ordering","aggregation: ordering")
      
         ("aggregation: drop scheme","aggregation: drop scheme")
      
         ("aggregation: threshold","aggregation: drop tol")
      
         ("aggregation: min agg size","aggregation: min agg size")
      
         ("aggregation: max agg size","aggregation: max agg size")
      
         ("aggregation: max selected neighbors","aggregation: max selected neighbors")
      
         ("aggregation: Dirichlet threshold","aggregation: Dirichlet threshold")
      
         ("aggregation: enable phase 1","aggregation: enable phase 1")
      
         ("aggregation: enable phase 2a","aggregation: enable phase 2a")
      
         ("aggregation: enable phase 2b","aggregation: enable phase 2b")
      
         ("aggregation: enable phase 3","aggregation: enable phase 3")
      
         ("aggregation: preserve Dirichlet points","aggregation: preserve Dirichlet points")
      
         ("aggregation: export visualization data","aggregation: export visualization data")
      
         ("aggregation: mode","aggregation: mode")
      
         ("export data","export data")
      
         ("ML print initial list","print initial parameters")
      
         ("print unused","print unused parameters")
      
         ("transpose: use implicit","transpose: use implicit")
      
         ("energy minimization: enable","multigrid algorithm")
      
         ("semicoarsen: coarsen rate","semicoarsen: coarsen rate")
      
         ("aggregation: damping factor","sa: damping factor")
      
         ("sa: use filtered matrix","sa: use filtered matrix")
      
         ("sa: calculate eigenvalue estimate","sa: calculate eigenvalue estimate")
      
         ("eigen-analysis: iterations","sa: eigenvalue estimate num iterations")
      
         ("filtered matrix: use lumping","filtered matrix: use lumping")
      
         ("filtered matrix: reuse eigenvalue","filtered matrix: reuse eigenvalue")
      
         ("filtered matrix: reuse graph","filtered matrix: reuse graph")
      
         ("emin: iterative method","emin: iterative method")
      
         ("emin: num iterations","emin: num iterations")
      
         ("emin: num reuse iterations","emin: num reuse iterations")
      
         ("emin: pattern","emin: pattern")
      
         ("emin: pattern order","emin: pattern order")
      
         ("repartition: enable","repartition: enable")
      
         ("repartition: partitioner","repartition: partitioner")
      
         ("repartition: params","repartition: params")
      
         ("repartition: start level","repartition: start level")
      
         ("repartition: min per proc","repartition: min rows per proc")
      
         ("repartition: max min ratio","repartition: max imbalance")
      
         ("repartition: remap parts","repartition: remap parts")
      
         ("repartition: remap num values","repartition: remap num values")
      
         ("repartition: print partition distribution","repartition: print partition distribution")
      
         ("repartition: rebalance P and R","repartition: rebalance P and R")
      
         ("reuse: type","reuse: type")
      ;

}

