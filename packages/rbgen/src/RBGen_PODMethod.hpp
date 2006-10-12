

#ifndef RBGEN_POD_METHOD_HPP
#define RBGEN_POD_METHOD_HPP

#include <vector>

namespace RBGen {

  template<class ScalarType>  
  class PODMethod {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    PODMethod() {};

    //! Destructor.
    virtual ~PODMethod() {};
    //@}

    //@{ @name Get methods
    
    //! Returns the singular values computed corresponding to the reduced basis.
    virtual std::vector<ScalarType> getSingularValues() const = 0;

    //@}
    
  };
  
} // end of RBGen namespace

#endif // RBGEN_POD_METHOD_HPP
