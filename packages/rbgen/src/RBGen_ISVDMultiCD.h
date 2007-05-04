#ifndef RBGEN_ISVDMULTICD_H
#define RBGEN_ISVDMULTICD_H

#include "RBGen_IncSVDPOD.h"

//
// computeBasis()
//       |
//   makePass()     ___
//   while () {    /    expand()
//      incStep()  ---- SVD()
//   }             \___ shrink()
//
// makePass(), expand() and shrink() are pure virtual
//
// makePass() is implemented in a base class that decides 
// if a method is Multipass, and if so, in what manner
// makePass() puts the data for the next pass into the proper columns
// of U_, then calls incstep (telling it how many new columns)
//
// incstep calls expand to construct expanded U_ and V_, and B_
// it computes the SVD of B_
// it passes this data to shrink(), which shrink according to the 
// base class.
//
// expand(),shrink() are implemented in a base class that 
// decides the representation: UDV, QRW, QBW
// 
//                                   IncSVDPOD
//                                       |
//      -------------------------------------------------------------------
//      |       |        |           |            |            |          |
//  ISVDUDV  ISVDQRW  ISVDQBW   ISVDMultiCD ISVDMultiSDA ISVDMultiSDB ISVDSingle
//      |       |        |           |            |            |          |
//      ------------------           --------------------------------------
//              \                                       /
//               \                                     /
//                \                                   /
//                 \                                 /
//                  \                               /
//                   \---  Concrete Base Class ----/
//
// Then a concrete base class (one of 3x4==12 varieties) is formed simply through
// inheritence. This is the Template Pattern type of Design Pattern.
//

namespace RBGen {

  //! Class for producing a basis using the Incremental SVD in a single pass
  class ISVDMultiCD : public virtual IncSVDPOD {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    ISVDMultiCD();

    //! Destructor.
    virtual ~ISVDMultiCD() {};
    //@}

    //! @name Set Methods
    //@{

    //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RefCountPtr< Teuchos::ParameterList >& params,
                     const Teuchos::RefCountPtr< Epetra_MultiVector >& init,
                     const Teuchos::RefCountPtr< RBGen::FileIOHandler< Epetra_CrsMatrix > >& fileio = Teuchos::null );

    //@}

  protected:

    // private member for performing inc steps
    virtual void makePass();

    // will need workspace for A*W = A - A Z T Z^T
    // * local multivector for Z
    // * dist multivector for A*Z
    Teuchos::RefCountPtr<Epetra_MultiVector> workAZT_, workZ_;
    Teuchos::RefCountPtr<Epetra_SerialDenseMatrix> workT_;

  };

} // end of RBGen namespace

#endif // RBGEN_ISVDMULTICD_H
