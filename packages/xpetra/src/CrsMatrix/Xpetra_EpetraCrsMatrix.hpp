#ifndef XPETRA_EPETRACRSMATRIX_HPP
#define XPETRA_EPETRACRSMATRIX_HPP

/* this file is automatically generated - do not edit (see script/epetra.py) */

#include "Xpetra_EpetraConfigDefs.hpp"

#include "Xpetra_CrsMatrix.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>

#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_EpetraVector.hpp"
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraCrsGraph.hpp"

#include "Xpetra_Utils.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  class EpetraCrsMatrix
    : public CrsMatrix<double, int, int>
  {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<void,int,Node>::SparseOps LocalMatOps;

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor specifying fixed number of entries for each row.
    EpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null);

    //! Constructor specifying (possibly different) number of entries in each row.
    EpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null);

    //! Constructor specifying column Map and fixed number of entries for each row.
    EpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null);

    //! Constructor specifying column Map and number of entries in each row.
    EpetraCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, ProfileType pftype=DynamicProfile, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null);

    //! Constructor specifying a previously constructed graph.
    EpetraCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > &graph, const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null);

    //! 
    virtual ~EpetraCrsMatrix() { }

    //@}

    //! @name Insertion/Removal Methods
    //@{

    //! Insert matrix entries, using global IDs.
    void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals);

    //! Scale the current values of a matrix, this = alpha*this.
    void scale(const Scalar &alpha) { mtx_->Scale(alpha); }

    //@}

    //! @name Transformational Methods
    //@{

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, OptimizeOption os=DoOptimizeStorage) { mtx_->FillComplete(toEpetra(domainMap), toEpetra(rangeMap), toEpetra(os)); }

    //! Signal that data entry is complete.
    void fillComplete(OptimizeOption os=DoOptimizeStorage) { mtx_->FillComplete(toEpetra(os)); }

    //@}

    //! @name Methods implementing RowMatrix
    //@{

    //! Returns the communicator.
    const RCP< const Comm< int > >  getComm() const { return toXpetra(mtx_->Comm()); }

    //! Returns the Map that describes the row distribution in this matrix.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const { return toXpetra(mtx_->RowMap()); }

    //! Returns the Map that describes the column distribution in this matrix.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const { return toXpetra(mtx_->ColMap()); }

    //! Returns the CrsGraph associated with this matrix.
    RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node, LocalMatOps > > getCrsGraph() const { return toXpetra(mtx_->Graph()); }

    //! Number of global elements in the row map of this matrix.
    global_size_t getGlobalNumRows() const { return mtx_->NumGlobalRows(); }

    //! Number of global columns in the matrix.
    global_size_t getGlobalNumCols() const { return mtx_->NumGlobalCols(); }

    //! Returns the number of matrix rows owned on the calling node.
    size_t getNodeNumRows() const { return mtx_->NumMyRows(); }

    //! Returns the number of columns connected to the locally owned rows of this matrix.
    size_t getNodeNumCols() const { return mtx_->NumMyCols(); }

    //! Returns the global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const { return mtx_->NumGlobalNonzeros(); }

    //! Returns the local number of entries in this matrix.
    size_t getNodeNumEntries() const { return mtx_->NumMyNonzeros(); }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const { return mtx_->NumMyEntries(localRow); }

    //! Returns the number of global diagonal entries, based on global row/column index comparisons.
    global_size_t getGlobalNumDiags() const { return mtx_->NumGlobalDiagonals(); }

    //! Returns the number of local diagonal entries, based on global row/column index comparisons.
    size_t getNodeNumDiags() const { return mtx_->NumMyDiagonals(); }

    //! Returns the maximum number of entries across all rows/columns on all nodes.
    size_t getGlobalMaxNumRowEntries() const { return mtx_->GlobalMaxNumEntries(); }

    //! Returns the maximum number of entries across all rows/columns on this node.
    size_t getNodeMaxNumRowEntries() const { return mtx_->MaxNumEntries(); }

    //! If matrix indices are in the local range, this function returns true. Otherwise, this function returns false.
    bool isLocallyIndexed() const { return mtx_->IndicesAreLocal(); }

    //! If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
    bool isGloballyIndexed() const { return mtx_->IndicesAreGlobal(); }

    //! Returns true if fillComplete() has been called and the matrix is in compute mode.
    bool isFillComplete() const { return mtx_->Filled(); }

    //! Returns the Frobenius norm of the matrix.
    ScalarTraits< Scalar >::magnitudeType getFrobeniusNorm() const { return mtx_->NormFrobenius(); }

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    void getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView< LocalOrdinal > &Indices, const ArrayView< Scalar > &Values, size_t &NumEntries) const;

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &indices, ArrayView< const Scalar > &values) const;

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices, ArrayView< const Scalar > &values) const;

    //! Get a copy of the diagonal entries owned by this node, with local row indices.
    void getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const { mtx_->ExtractDiagonalCopy(toEpetra(diag)); }

    //@}

    //! @name Methods implementing Operator
    //@{

    //! Computes the sparse matrix-multivector multiplication.
    void apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=ScalarTraits< Scalar >::one(), Scalar beta=ScalarTraits< Scalar >::zero()) const;

    //! Returns the Map associated with the domain of this operator. This will be null until fillComplete() is called.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const { return toXpetra(mtx_->DomainMap()); }

    //! 
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const { return toXpetra(mtx_->RangeMap()); }

    //@}

    //! @name Overridden from Teuchos::Describable
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const { return toXpetra(mtx_->Map()); }

    //! Import.
    void doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source, const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM);

    //! Export.
    void doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest, const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM);

    //! Import (using an Exporter).
    void doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source, const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM);

    //! Export (using an Importer).
    void doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest, const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM);

    //@}

    //! @name Xpetra specific
    //@{

    //! EpetraCrsMatrix constructor to wrap a Epetra_CrsMatrix object
    EpetraCrsMatrix(const Teuchos::RCP<Epetra_CrsMatrix > &mtx) : mtx_(mtx) {  }

    //! Get the underlying Epetra matrix
    RCP<const Epetra_CrsMatrix> getEpetra_CrsMatrix() const { return mtx_; }
    
    //! Get the underlying Epetra matrix
    RCP<Epetra_CrsMatrix> getEpetra_CrsMatrixNonConst() const { return mtx_; } //TODO: remove
 
   //@}
    
  private:
    
    RCP<Epetra_CrsMatrix> mtx_;

  }; // EpetraImport class

} // Xpetra namespace

#define XPETRA_EPETRACRSMATRIX_SHORT
#endif // XPETRA_EPETRACRSMATRIX_HPP
