// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
// @HEADER

#ifndef TPETRA_CRSMATRIX_DECL_HPP
#define TPETRA_CRSMATRIX_DECL_HPP

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient

// TODO: add typeglobs: CrsMatrix<Scalar,typeglob>
// TODO: add template (template) parameter for nonlocal container (this will be part of typeglob)

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsMatrix.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_decl.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
  // struct for i,j,v triplets
  template <class Ordinal, class Scalar>
  struct CrsIJV {
    CrsIJV();
    CrsIJV(Ordinal row, Ordinal col, const Scalar &val);
    Ordinal i,j;
    Scalar  v;
  };
}

namespace Teuchos {
  // SerializationTraits specialization for CrsIJV, using DirectSerialization
  template <typename Ordinal, typename Scalar>
  class SerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  : public DirectSerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  {};
}

namespace std {
  template <class Ordinal, class Scalar>
  bool operator<(const Tpetra::CrsIJV<Ordinal,Scalar> &ijv1, const Tpetra::CrsIJV<Ordinal,Scalar> &ijv2);
}
#endif

namespace Tpetra {

  //! \brief Sparse matrix that presents a compressed sparse row interface.
  /*!
   \tparam Scalar The type of the numerical entries of the matrix.
     (You can use real-valued or complex-valued types here, unlike in
     Epetra, where the scalar type is always \c double.)

   \tparam LocalOrdinal The type of local indices.  Same as the \c
     LocalOrdinal template parameter of \c Map objects used by this
     matrix.  (In Epetra, this is just \c int.)

   \tparam GlobalOrdinal The type of global indices.  Same as the \c
     GlobalOrdinal template parameter of \c Map objects used by this
     matrix.  (In Epetra, this is just \c int.  One advantage of
     Tpetra over Epetra is that you can use a 64-bit integer type here
     if you want to solve big problems.)

   \tparam Node A class implementing on-node shared-memory parallel
     operations.  It must implement the \ref kokkos_node_api "Kokkos
     Node API."  The default \c Node type depends on your build
     options.

   \tparam LocalMatOps A local sparse matrix operations class.  It
     must implement the \ref kokkos_crs_ops "Kokkos CRS Ops API."

   This class implements a distributed-memory parallel sparse matrix,
   and provides sparse matrix-vector multiply (including transpose)
   and sparse triangular solve operations.  It provides access by rows
   to the elements of the matrix, as if the local data were stored in
   compressed sparse row format.  This class has an interface like
   that of \c Epetra_CrsMatrix, but also allows insertion of data into
   nonowned rows, much like Epetra_FECrsMatrix.

   <b>Local vs. Global</b>

   The distinction between local and global indices might confuse new
   Tpetra users.  <it>Global</it> indices represent the rows and
   columns uniquely over the entire matrix, which may be distributed
   over multiple processes.  <it>Local</it> indices are local to the
   process that owns them.  If global index G is owned by process P,
   then there is a unique local index L on process P corresponding to
   G.  If the local index L is valid on process P, then there is a
   unique global index G owned by P corresponding to the pair (L, P).
   However, multiple processes might own the same global index (an
   "overlapping Map"), so a global index G might correspond to
   multiple (L, P) pairs.  In summary, local indices on a process
   correspond to matrix rows or columns owned by that process.

   We summarize the different between local and global indices because
   many of CrsMatrix's methods for adding, modifying, or accessing
   entries come in versions that take either local or global indices.
   The matrix itself may store indices either as local or global.  You
   should only use the method version corresponding to the current
   state of the matrix.  For example, \c getGlobalRowView() returns a
   view to the indices represented as global; it is incorrect to call
   this method if the matrix is storing indices as local.  Call the \c
   isGloballyIndexed() or \c isLocallyIndexed() methods to find out
   whether the matrix currently stores indices as local or global.

   Method methods that work with global indices only allow operations
   on indices owned by the calling process.  For example, methods that
   take a global row index expect that row to be owned by the calling
   process.  Access to nonlocal (i.e., not owned by the calling
   process) rows requires performing an explicit communication via the
   import/export capabilities of the CrsMatrix object; see \c
   DistObject.  However, the method \c insertGlobalValues() is an
   exception to this rule.  It allows you to add data to nonlocal
   rows.  These data are stored locally and communicated to the
   appropriate node on the next call to \c globalAssemble() or \c
   fillComplete() (the latter calls the former).
  */
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >
  class CrsMatrix : public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
                    public DistObject<char, LocalOrdinal,GlobalOrdinal,Node> {
    public:
      typedef Scalar        scalar_type;
      typedef LocalOrdinal  local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef Node          node_type;
      // backwards compatibility defines both of these
      typedef LocalMatOps   mat_vec_type;
      typedef LocalMatOps   mat_solve_type;

      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor specifying the number of non-zeros for all rows.
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor specifying the number of non-zeros for each row.
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a column map and the number of non-zeros for all rows.
      /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
        */
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a column map and the number of non-zeros for each row.
      /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
        */
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a pre-constructed graph.
      explicit CrsMatrix(const RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &graph);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! Insert matrix entries, using global IDs.
      /** All index values must be in the global space. 
          \pre \c globalRow exists as an ID in the global row map
          \pre <tt>isStorageOptimized() == false</tt>

          \note If \c globalRow does not belong to the matrix on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local matrix. 
          \note If the matrix row already contains values at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
          \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
          \note If <tt>isLocallyIndexed() == true</tt>, then the global indices will be translated to local indices via the column map; indices not present in the column map will be discarded.
       */
      void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals);

      //! Insert matrix entries, using local IDs.
      /**
          \pre \c localRow is a local row belonging to the matrix on this node
          \pre <tt>isGloballyIndexed() == false</tt>
          \pre <tt>isStorageOptimized() == false</tt>
          \pre <tt>hasColMap() == true</tt>

          \post <tt>isLocallyIndexed() == true</tt>

          \note If the matrix row already contains entries at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
          \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
        */
      void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals);

      //! \brief Replace matrix entries, using global IDs.
      /** All index values must be in the global space. 

         \pre \c globalRow is a global row belonging to the matrix on this node.

         \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
      void replaceGlobalValues(GlobalOrdinal globalRow, 
                               const ArrayView<const GlobalOrdinal> &cols,
                               const ArrayView<const Scalar>        &vals);

      //! Replace matrix entries, using local IDs.
      /** All index values must be in the local space. 
        */
      void replaceLocalValues(LocalOrdinal localRow, 
                              const ArrayView<const LocalOrdinal> &cols,
                              const ArrayView<const Scalar>       &vals);

      //! Sum into multiple entries, using global IDs.
      /** All index values must be in the global space. 

         \pre \c globalRow is a global row belonging to the matrix on this node.

        */
      void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                               const ArrayView<const GlobalOrdinal> &cols,
                               const ArrayView<const Scalar>        &vals);


      //! Sum into multiple entries, using local IDs.
      /** All index values must be in the local space. 

         \pre \c localRow is a local row belonging to the matrix on this node.

        */
      void sumIntoLocalValues(LocalOrdinal globalRow, 
                              const ArrayView<const LocalOrdinal>  &cols,
                              const ArrayView<const Scalar>        &vals); 

      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(const Scalar &alpha);

      //! Scale the current values of a matrix, this = alpha*this. 
      void scale(const Scalar &alpha);

      //@}
	  

      //! @name Transformational Methods
      //@{ 

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      /*! Resume fill operations.
          After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

          resumeFill() may be called repeatedly. 

          \post  <tt>isFillActive() == true<tt>
          \post  <tt>isFillComplete() == false<tt>
       */
      void resumeFill();

      /*! \brief Signal that data entry is complete, specifying domain and range maps.

          Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

          \pre  <tt>isFillActive() == true<tt>
          \pre <tt>isFillComplete()() == false<tt>

          \post <tt>isFillActive() == false<tt>
          \post <tt>isFillComplete() == true<tt>
          \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>. See isStorageOptimized() for consequences.
       */ 
      void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage);

      /*! \brief Signal that data entry is complete. 

          Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

          \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

          \pre  <tt>isFillActive() == true<tt>
          \pre <tt>isFillComplete()() == false<tt>

          \post <tt>isFillActive() == false<tt>
          \post <tt>isFillComplete() == true<tt>
          \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>. See isStorageOptimized() for consequences.
       */
      void fillComplete(OptimizeOption os = DoOptimizeStorage);

      //@}

      //! @name Methods implementing RowMatrix
      //@{ 

      //! Returns the communicator.
      const RCP<const Comm<int> > & getComm() const;

      //! Returns the underlying node.
      RCP<Node> getNode() const;

      //! Returns the Map that describes the row distribution in this matrix.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

      //! Returns the RowGraph associated with this matrix. 
      RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const;

      //! Returns the CrsGraph associated with this matrix. 
      RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > getCrsGraph() const;

      /// \brief Number of global elements in the row map of this matrix.
      ///
      /// This is <it>not</it> the number of rows in the matrix as a
      /// mathematical object.  This method returns the global sum of
      /// the number of local elements in the row map on each
      /// processor, which is the row map's getGlobalNumElements().
      /// Since the row map is not one-to-one in general, that global
      /// sum could be different than the number of rows in the
      /// matrix.  If you want the number of rows in the matrix, ask
      /// the range map for its global number of elements, using the
      /// following code:
      /// <code>
      /// global_size_t globalNumRows = getRangeMap()->getGlobalNumElements();
      /// </code>
      /// This method retains the behavior of Epetra, which also asks
      /// the row map for the global number of rows, rather than
      /// asking the range map.
      ///
      /// \warning Undefined if isFillActive().
      /// 
      global_size_t getGlobalNumRows() const;

      //! \brief Number of global columns in the matrix.
      /** Returns the number of entries in the domain map of the matrix.
          \warning Undefined if isFillActive().
      */
      global_size_t getGlobalNumCols() const;

      //! Returns the number of matrix rows owned on the calling node.
      size_t getNodeNumRows() const;

      //! Returns the number of columns connected to the locally owned rows of this matrix.
      /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
        */
      size_t getNodeNumCols() const;

      //! Returns the index base for global indices for this matrix. 
      GlobalOrdinal getIndexBase() const;

      //! Returns the global number of entries in this matrix.
      global_size_t getGlobalNumEntries() const;

      //! Returns the local number of entries in this matrix.
      size_t getNodeNumEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this matrix. */
      size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
      size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      /** Undefined if isFillActive().
        */
      global_size_t getGlobalNumDiags() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      /** Undefined if isFillActive().
        */
      size_t getNodeNumDiags() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
      /** Undefined if isFillActive().
        */
      size_t getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node.
      /** Undefined if isFillActive().
        */
      size_t getNodeMaxNumRowEntries() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      bool hasColMap() const; 

      //! \brief Indicates whether the matrix is lower triangular.
      /** Undefined if isFillActive().
        */
      bool isLowerTriangular() const;

      //! \brief Indicates whether the matrix is upper triangular.
      /** Undefined if isFillActive().
        */
      bool isUpperTriangular() const;

      //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false.
      bool isLocallyIndexed() const;

      //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
      bool isGloballyIndexed() const;

      //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
      bool isFillComplete() const;

      //! Returns \c true if resumeFill() has been called and the matrix is in edit mode.
      bool isFillActive() const;

      //! \brief Returns \c true if storage has been optimized.
      /**
        Optimized storage means that the allocation of each row is equal to the
        number of entries. The effect is that a pass through the matrix, i.e.,
        during a mat-vec, requires minimal memory traffic. One limitation of
        optimized storage is that no new indices can be added to the matrix.
        */
      bool isStorageOptimized() const;

      //! Returns \c true if the matrix was allocated with static data structures.
      ProfileType getProfileType() const;

      //! Indicates that the graph is static, so that new entries cannot be added to this matrix.
      bool isStaticGraph() const;

      //! Returns the Frobenius norm of the matrix. 
      /** Computes and returns the Frobenius norm of the matrix, defined as:
        \f$ \|A\|_F = \sqrt{\sum_{i,j} \|\a_{ij}\|^2} \f$

         If the matrix is fill-complete, then the computed value is cached; the cache is cleared whenever resumeFill() is called.
         Otherwise, the value is computed every time the method is called.
        */
      typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const;

      //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param Values - (Out) Matrix values.
        \param NumEntries - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
         returned as OrdinalTraits<size_t>::invalid().
       */
      void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                            const ArrayView<GlobalOrdinal> &Indices,
                            const ArrayView<Scalar> &Values,
                            size_t &NumEntries
                            ) const;

      //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices - (Out) Local column indices corresponding to values.
        \param Values - (Out) Matrix values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
         with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
         returned as OrdinalTraits<size_t>::invalid().

        \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
       */
      void getLocalRowCopy(LocalOrdinal LocalRow, 
                           const ArrayView<LocalOrdinal> &Indices, 
                           const ArrayView<Scalar> &Values,
                           size_t &NumEntries
                           ) const;

      //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
      /*!
        \param GlobalRow - (In) Global row number for which indices are desired.
        \param Indices   - (Out) Global column indices corresponding to values.
        \param Values    - (Out) Row values
        \pre <tt>isLocallyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

         Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
       */
      void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const;

      //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices  - (Out) Global column indices corresponding to values.
        \param Values   - (Out) Row values
        \pre <tt>isGloballyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

         Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
       */
      void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const;

      //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
      /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const;

      /** \brief . */
      void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);

      /** \brief . */
      void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);

      //@}

      //! @name Advanced templated methods
      //@{

      //! Multiplies this matrix by a MultiVector.
      /*! If \c trans is \c Teuchos::NO_TRANS, then  X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix. 
          Otherwise, then  X is should be described by the row map of the matrix and \c Y should be described by the column map of the matrix. 

          Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

          This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
          If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
          will be accumulated into \c Y.
       */
      template <class DomainScalar, class RangeScalar>
      void localMultiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;

      //! Solves a linear system when the underlying matrix is triangular.
      /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

          This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that solve() not be called directly; instead, use the CrsMatrixSolveOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
          Both are required to have constant stride. However, unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No runtime checking will be performed in a non-debug build.
       */
      template <class DomainScalar, class RangeScalar>
      void localSolve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;
          
      //@}

      //! Returns another CrsMatrix with the same entries, but represented as a different scalar type.
      template <class T>
      RCP<CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node> > convert() const;
          
      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! \brief Computes the sparse matrix-multivector multiplication.
      /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
          - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
       */
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = ScalarTraits<Scalar>::one(),
                 Scalar beta = ScalarTraits<Scalar>::zero()) const;

      //! Indicates whether this operator supports applying the adjoint operator.
      bool hasTransposeApply() const;

      //! \brief Returns the Map associated with the domain of this operator.
      //! This will be <tt>null</tt> until fillComplete() is called.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This will be <tt>null</tt> until fillComplete() is called.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

      //@}

      //! @name Overridden from Teuchos::Describable 
      //@{

      /** \brief Return a simple one-line description of this object. */
      std::string description() const;

      /** \brief Print the object with some verbosity level to an FancyOStream object. */
      void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

      //@}

      //! @name Methods implementing Tpetra::DistObject
      //@{

      bool checkSizes(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source);

      void copyAndPermute(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                          size_t numSameIDs,
                          const ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const ArrayView<const LocalOrdinal> &permuteFromLIDs);

      void packAndPrepare(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                          const ArrayView<const LocalOrdinal> &exportLIDs,
                          Array<char> &exports,
                          const ArrayView<size_t> & numPacketsPerLID,
                          size_t& constantNumPackets,
                          Distributor &distor);

      /// \brief Unpack the imported column indices and values, and combine into matrix.
      ///
      /// \warning The allowed \c combineMode depends on whether the
      ///   matrix's graph is static or dynamic.  ADD, REPLACE, and
      ///   ABSMAX are valid for a static graph, but INSERT is not.
      ///   ADD and INSERT are valid for a dynamic graph; ABSMAX and
      ///   REPLACE have not yet been implemented (and would require
      ///   serious changes to matrix assembly in order to implement
      ///   sensibly).
      void 
      unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
			const Teuchos::ArrayView<const char> &imports,
			const Teuchos::ArrayView<size_t> &numPacketsPerLID,
			size_t constantNumPackets,
			Distributor &distor,
			CombineMode combineMode);
      //@}

      //! \name Deprecated routines to be removed at some point in the future.
      //@{

      /** \brief Deprecated. Re-allocate the data into contiguous storage.

          This method is deprecated and will be removed in a future version of Tpetra, as 
          the implementation of storage optimization has been below Tpetra to Kokkos.

          Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
          required to be called by all nodes that participate in the associated communicator.
       */
      TPETRA_DEPRECATED void optimizeStorage();

      //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
      TPETRA_DEPRECATED void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayRCP<const GlobalOrdinal> &indices, ArrayRCP<const Scalar> &values) const;

      //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
      TPETRA_DEPRECATED void getLocalRowView(LocalOrdinal LocalRow, ArrayRCP<const LocalOrdinal> &indices, ArrayRCP<const Scalar> &values) const;

      //! Deprecated. Replaced by localMultiply().
      template <class DomainScalar, class RangeScalar>
      TPETRA_DEPRECATED 
      void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;

      //! Deprecated. Replaced by localSolve().
      template <class DomainScalar, class RangeScalar>
      TPETRA_DEPRECATED 
      void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;

      //@}

    private:
      // We forbid copy construction by declaring this method private
      // and not implementing it.
      CrsMatrix (const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &rhs);

      // We forbid assignment (operator=) by declaring this method
      // private and not implementing it.
      CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>& 
      operator= (const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &rhs);

      /// \brief Combine in the data using the given combine mode.
      ///
      /// The \c copyAndPermute() and \c unpackAndCombine() methods
      /// use this function to combine incoming entries from the
      /// source matrix with the target matrix's current data.  This
      /// method's behavior depends on whether the target matrix (that
      /// is, this matrix) has a static graph.
      void 
      combineGlobalValues (const GlobalOrdinal globalRowIndex, 
			   const Teuchos::ArrayView<const GlobalOrdinal> columnIndices,
			   const Teuchos::ArrayView<const Scalar> values,
			   const Tpetra::CombineMode combineMode);

      /// \brief Transform CrsMatrix entries by applying a binary function to them.
      ///
      /// For every entry \f$A(i,j)\f$ to transform, if \f$v_{ij}\f$ is
      /// the corresponding entry of the \c values array, then we apply
      /// the function to \f$A(i,j)\f$ as follows:
      /// \f[
      ///   A(i,j) := f(A(i,j), v_{ij}).
      /// \f]
      /// For example, BinaryFunction = std::plus<Scalar> implements
      /// \c sumIntoGlobalValues(), and BinaryFunction =
      /// secondArg<Scalar,Scalar> implements replaceGlobalValues().
      ///
      /// \tparam BinaryFunction The type of binary function to apply.
      ///   std::binary_function is a model for this.
      template<class BinaryFunction>
      void
      transformGlobalValues (GlobalOrdinal globalRow, 
			     const Teuchos::ArrayView<const GlobalOrdinal>& indices,
			     const Teuchos::ArrayView<const Scalar>        & values,
			     BinaryFunction f)
      {
	typedef Scalar ST;
	typedef LocalOrdinal LO;
	typedef GlobalOrdinal GO;
	typedef Node NT;
	using Teuchos::Array;
	using Teuchos::ArrayView;

	TEUCHOS_TEST_FOR_EXCEPTION(values.size() != indices.size(), 
          std::invalid_argument, "transformGlobalValues: values.size() = " 
	  << values.size() << " != indices.size() = " << indices.size() << ".");

	const LO lrow = this->getRowMap()->getLocalElement(globalRow);

	TEUCHOS_TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::invalid_argument, 
	  "transformGlobalValues: The given global row index " << globalRow 
	  << " is not owned by the calling process (rank " 
	  << this->getRowMap()->getComm()->getRank() << ").");

	RowInfo rowInfo = staticGraph_->getRowInfo(lrow);
	if (indices.size() > 0) {
	  if (isLocallyIndexed()) {
	    // Convert global indices to local indices.
	    const Map<LO, GO, NT> &colMap = *(this->getColMap());
	    Array<LO> lindices (indices.size());
	    typename ArrayView<const GO>::iterator gindit = indices.begin();
	    typename Array<LO>::iterator           lindit = lindices.begin();
	    while (gindit != indices.end()) {
	      // No need to filter before asking the column Map to
	      // convert GID->LID.  If the GID doesn't exist in the
	      // column Map, the GID will be mapped to invalid(), which
	      // will not be found in the graph.
	      *lindit++ = colMap.getLocalElement(*gindit++);
	    }
	    typename Graph::SLocalGlobalViews inds_view;
	    inds_view.linds = lindices();
	    staticGraph_->template transformValues<LocalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), f);
	  }
	  else if (isGloballyIndexed()) {
	    typename Graph::SLocalGlobalViews inds_view;
	    inds_view.ginds = indices;
	    staticGraph_->template transformValues<GlobalIndices>(rowInfo, inds_view, this->getViewNonConst(rowInfo).begin(), values.begin(), f);
	  }
	}
      }

    protected:
      // useful typedefs
      typedef OrdinalTraits<LocalOrdinal>                     LOT;
      typedef OrdinalTraits<GlobalOrdinal>                    GOT;
      typedef ScalarTraits<Scalar>                             ST;
      typedef typename ST::magnitudeType                Magnitude;
      typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>       MV;
      typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>             V;
      typedef CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>  Graph;
      // Enums
      enum GraphAllocationStatus {
        GraphAlreadyAllocated,
        GraphNotYetAllocated
      };

      /// \brief Allocate values (and optionally indices) using the Node.
      ///
      /// \param gas [in] If GraphNotYetAllocated, allocate the
      ///   indices of \c myGraph_ via \c allocateIndices(lg) before
      ///   allocating values.
      ///
      /// \param lg [in] Argument passed into \c
      ///   myGraph_->allocateIndices(), if applicable.
      ///
      /// \pre If the graph (that is, staticGraph_) indices are
      ///   already allocated, then gas must be GraphAlreadyAllocated.
      ///   Otherwise, gas must be GraphNotYetAllocated.  We only
      ///   check for this precondition in debug mode.
      ///
      /// \pre If the graph indices are not already allocated, then
      ///   the graph must be owned by the matrix.
      void allocateValues (ELocalGlobal lg, GraphAllocationStatus gas);

      // Sorting and merging
      void sortEntries();
      void mergeRedundantEntries();
      // global consts
      void clearGlobalConstants();
      void computeGlobalConstants();
      // matrix data accessors
      ArrayView<const Scalar>    getView(RowInfo rowinfo) const;
      ArrayView<      Scalar>    getViewNonConst(RowInfo rowinfo);
      // local Kokkos objects
      void pushToLocalMatrix();
      void pullFromLocalMatrix();
      void fillLocalMatrix(OptimizeOption os);
      void fillLocalSparseOps();
      // debugging
      void checkInternalState() const;

      // Two graph pointers needed in order to maintain const-correctness:
      // staticGraph_ is a graph passed to the constructor. We are not allowed to modify it. it is always a valid pointer.
      // myGraph_     is a graph created here. We are allowed to modify it. if myGraph_ != null, then staticGraph_ = myGraph_
      RCP<const Graph> staticGraph_;
      RCP<      Graph>     myGraph_;

      Kokkos::CrsMatrix<Scalar,LocalOrdinal,Node,LocalMatOps> lclMatrix_;
      typename LocalMatOps::template rebind<Scalar>::other    lclMatOps_;

      // matrix values. before allocation, both are null.
      // after allocation, one is null.
      // 1D == StaticAllocation, 2D == DynamicAllocation
      // The allocation always matches that of graph_, as the graph does the allocation for the matrix.
      ArrayRCP<Scalar>                       values1D_;
      ArrayRCP<ArrayRCP<Scalar> >            values2D_;
      // TODO: these could be allocated at resumeFill() and de-allocated at fillComplete() to make for very fast getView()/getViewNonConst()
      // ArrayRCP< typedef ArrayRCP<const Scalar>::iterator > rowPtrs_;
      // ArrayRCP< typedef ArrayRCP<      Scalar>::iterator > rowPtrsNC_;

      bool fillComplete_;

      // non-local data
      std::map<GlobalOrdinal, Array<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

      // a wrapper around multiply, for use in apply; it contains a non-owning RCP to *this, therefore, it is not allowed 
      // to persist past the destruction of *this. therefore, WE MAY NOT SHARE THIS POINTER.
      RCP< const CrsMatrixMultiplyOp<Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > sameScalarMultiplyOp_;

      // cached frobenius norm: -ST::one() means invalid
      mutable Magnitude frobNorm_;

  }; // class CrsMatrix

  /** \brief Non-member function to create an empty CrsMatrix given a row map and a non-zero profile.

      Returns a dynamically allocated (DynamicProfile) matrix with specified number of non-zeros per row (defaults to zero).

      \relatesalso CrsMatrix
   */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createCrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t maxNumEntriesPerRow = 0)
  {
    RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ret;
    ret = rcp(new CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map, maxNumEntriesPerRow, DynamicProfile) );
    return ret;
  }

} // namespace Tpetra

/**
  \example LocalMatOpExample.cpp
  An example using a different sparse mat-vec with Tpetra::CrsMatrix and Tpetra::CrsGraph.
 */

/** 
  \example CrsMatrix_NonlocalAfterResume.hpp
  An example for inserting non-local entries into a Tpetra::CrsMatrix using Tpetra::CrsMatrix::insertGlobalValues(), with multiple calls to Tpetra::CrsMatrix::fillComplete().
 */

#endif
