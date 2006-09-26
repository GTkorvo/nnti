

#ifndef RBGEN_MS_PREPROCESSOR_H
#define RBGEN_MS_PREPROCESSOR_H

#include "RBGen_Preprocessor.hpp"
#include "RBGen_FileIOHandler.hpp"

// Forward declaration of Epetra_MultiVector class
class Epetra_MultiVector;

namespace RBGen {
  
  class MSPreprocessor : public Preprocessor< Epetra_MultiVector > {
    
  public:
    //@{ @name Constructor/Destructor.

    //! Default constructor.
    MSPreprocessor();

    //! Destructor.
    virtual ~MSPreprocessor() {};
    //@}

    //@{ @name Initialization/Reset Methods

    //! Initialize preprocessor
    void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params, 
                     const Teuchos::RefCountPtr< FileIOHandler <Epetra_MultiVector> >& fileio );

    void Reset() { isInitialized_ = false; };
    //@}

    //@{ @name Preprocess Methods

    //! Preprocess the snapshot set passed in
    void Preprocess( Teuchos::RefCountPtr<Epetra_MultiVector>& ss );
    //@}

    //@{ @name Status Methods

    //! Return initialized status of the preprocessor
    bool isInitialized() const { return isInitialized_; };

    //@}
  private:

    //! Initialization flag.
    bool isInitialized_;

    //! Steady state filename.
    std::string steady_file_;

    //! Scaling vector for the snapshots.
    std::vector< double > scalings_;
    
    //! Scaling indices for the snapshots.
    std::vector< std::pair<int,int> > scaling_idx_;

    //! Pointer to the File I/O Handler object.
    Teuchos::RefCountPtr< FileIOHandler< Epetra_MultiVector > > fileio_;
  };
  
} // end of RBGen namespace

#endif // RBGEN_MS_PREPROCESSOR_H
