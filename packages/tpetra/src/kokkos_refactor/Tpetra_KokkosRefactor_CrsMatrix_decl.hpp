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

#ifndef TPETRA_KOKKOSREFACTOR_CRSMATRIX_DECL_HPP
#define TPETRA_KOKKOSREFACTOR_CRSMATRIX_DECL_HPP

/*#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix_decl.hpp"
#include "Tpetra_Exceptions.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_KokkosRefactor_CrsGraph.hpp"
#include "Tpetra_KokkosRefactor_Vector.hpp"*/
//#include "Tpetra_CrsMatrix_decl.hpp"
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

#include "Kokkos_CrsMatrix.hpp"
namespace Tpetra {

  //! \brief Sparse matrix that presents a compressed sparse row interface.
  /*!
   \tparam Scalar The type of the numerical entries of the matrix.
     (You can use real-valued or complex-valued types here, unlike in
     Epetra, where the scalar type is always \c double.)

   \tparam LocalOrdinal The type of local indices.  Same as the
     <tt>LocalOrdinal</tt> template parameter of Map objects used by
     this matrix.  (In Epetra, this is just \c int.)  The default type
     is \c int, which should suffice for most users.  This type must
     be big enough to store the local (per process) number of rows or
     columns.

   \tparam GlobalOrdinal The type of global indices.  Same as the
     <tt>GlobalOrdinal</tt> template parameter of Map objects used by
     this matrix.  (In Epetra, this is just \c int.  One advantage of
     Tpetra over Epetra is that you can use a 64-bit integer type here
     if you want to solve big problems.)  The default type is
     <tt>LocalOrdinal</tt>.  This type must be big enough to store the
     global (over all processes in the communicator) number of rows or
     columns.

   \tparam Node A class implementing on-node shared-memory parallel
     operations.  It must implement the
     \ref kokkos_node_api "Kokkos Node API."
     The default \c Node type should suffice for most users.
     The actual default type depends on your Trilinos build options.

   \tparam LocalMatOps Type implementing local sparse
     matrix-(multi)vector multiply and local sparse triangular solve.
     It must implement the \ref kokkos_crs_ops "Kokkos CRS Ops API."
     The default \c LocalMatOps type should suffice for most users.
     The actual default type depends on your Trilinos build options.

   \note If you use the default \c GlobalOrdinal type, which is
     <tt>int</tt>, then the <i>global</i> number of rows or columns in
     the matrix may be no more than \c INT_MAX, which for typical
     32-bit \c int is \f$2^{31} - 1\f$ (about two billion).  If you
     want to solve larger problems, you must use a 64-bit integer type
     here.

   This class implements a distributed-memory parallel sparse matrix,
   and provides sparse matrix-vector multiply (including transpose)
   and sparse triangular solve operations.  It provides access by rows
   to the elements of the matrix, as if the local data were stored in
   compressed sparse row format.  (Implementations are <i>not</i>
   required to store the data in this way internally.)  This class has
   an interface like that of Epetra_CrsMatrix, but also allows
   insertion of data into nonowned rows, much like Epetra_FECrsMatrix.

   \section Tpetra_CrsMatrix_prereq Prerequisites

   Before reading the rest of this documentation, it helps to know
   something about the Teuchos memory management classes, in
   particular Teuchos::RCP, Teuchos::ArrayRCP, and Teuchos::ArrayView.
   You should also know a little bit about MPI (the Message Passing
   Interface for distributed-memory programming).  You won't have to
   use MPI directly to use CrsMatrix, but it helps to be familiar with
   the general idea of distributed storage of data over a
   communicator.  Finally, you should read the documentation of Map
   and MultiVector.

   \section Tpetra_CrsMatrix_local_vs_global Local and global indices

   The distinction between local and global indices might confuse new
   Tpetra users.  Please refer to the documentation of Map for a
   detailed explanation.  This is important because many of
   CrsMatrix's methods for adding, modifying, or accessing entries
   come in versions that take either local or global indices.  The
   matrix itself may store indices either as local or global, and the
   same matrix may use global indices or local indices at different
   points in its life.  You should only use the method version
   corresponding to the current state of the matrix.  For example,
   getGlobalRowView() returns a view to the indices represented as
   global; it is incorrect to call this method if the matrix is
   storing indices as local.  Call isGloballyIndexed() or
   isLocallyIndexed() to find out whether the matrix currently stores
   indices as local or global.

   \section Tpetra_CrsMatrix_insertion_into_nonowned_rows Insertion into nonowned rows

   All methods (except for insertGlobalValues() and
   sumIntoGlobalValues(); see below) that work with global indices
   only allow operations on indices owned by the calling process.  For
   example, methods that take a global row index expect that row to be
   owned by the calling process.  Access to <i>nonowned rows</i>, that
   is, rows <i>not</i> owned by the calling process, requires
   performing an explicit communication via the Import / Export
   capabilities of the CrsMatrix object.  See the documentation of
   DistObject for more details.

   The methods insertGlobalValues() and sumIntoGlobalValues() are
   exceptions to this rule.  They both allows you to add data to
   nonowned rows.  These data are stored locally and communicated to
   the appropriate process on the next call to globalAssemble() or
   fillComplete().  This means that CrsMatrix provides the same
   nonowned insertion functionality that Epetra provides via
   Epetra_FECrsMatrix.

   \section Tpetra_DistObject_MultDist Note for developers on DistObject

   DistObject only takes a single Map as input to its constructor.
   MultiVector is an example of a subclass for which a single Map
   suffices to describe its data distribution.  In that case,
   DistObject's getMap() method obviously must return that Map.
   CrsMatrix is an example of a subclass that requires two Map
   objects: a row Map and a column Map.  For CrsMatrix, getMap()
   returns the row Map.  This means that doTransfer() (which CrsMatrix
   does not override) uses the row Map objects of the source and
   target CrsMatrix objects.  CrsMatrix in turn uses its column Map
   (if it has one) to "filter" incoming sparse matrix entries whose
   column indices are not in that process' column Map.  This means
   that CrsMatrix may perform extra communication, though the Import
   and Export operations are still correct.

   This is necessary if the CrsMatrix does not yet have a column Map.
   Other processes might have added new entries to the matrix; the
   calling process has to see them in order to accept them.  However,
   the CrsMatrix may already have a column Map, for example, if it was
   created with the constructor that takes both a row and a column
   Map, or if it is fill complete (which creates the column Map if the
   matrix does not yet have one).  In this case, it could be possible
   to "filter" on the sender (instead of on the receiver, as CrsMatrix
   currently does) and avoid sending data corresponding to columns
   that the receiver does not own.  Doing this would require revising
   the Import or Export object (instead of the incoming data) using
   the column Map, to remove global indices and their target process
   ranks from the send lists if the target process does not own those
   columns, and to remove global indices and their source process
   ranks from the receive lists if the calling process does not own
   those columns.  (Abstractly, this is a kind of set difference
   between an Import or Export object for the row Maps, and the Import
   resp. Export object for the column Maps.)  This could be done
   separate from DistObject, by creating a new "filtered" Import or
   Export object, that keeps the same source and target Map objects
   but has a different communication plan.  We have not yet
   implemented this optimization.
  */
  template <class Scalar,
            class LocalOrdinal  ,
            class GlobalOrdinal ,
            class DeviceType>
  class CrsMatrix<Scalar, LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<Scalar,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps> :
                    public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >,
                    public DistObject<char, LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > {
  public:
    //! @name Typedefs
    //@{

    //! This class' first template parameter; the type of entries in the matrix.
    typedef Scalar                                scalar_type;
    //! This class' second template parameter; the type of local indices.
    typedef LocalOrdinal                          local_ordinal_type;
    //! This class' third template parameter; the type of global indices.
    typedef GlobalOrdinal                         global_ordinal_type;
    //! This class' fourth template parameter; the Kokkos Node type.
    typedef Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>  Node;
    typedef Node                                  node_type;
    /// \brief This class' fifth template parameter; the implementation of local sparse kernels.
    ///
    /// We define both this typedef and mat_solve_type for backwards compatibility.
    typedef typename KokkosClassic::DefaultKernels<Scalar,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps LocalMatOps;
    typedef LocalMatOps   mat_vec_type;
    /// \brief This class' fifth template parameter; the implementation of local sparse kernels.
    ///
    /// We define both this typedef and mat_vec_type for backwards compatibility.
    typedef LocalMatOps   mat_solve_type;

    //! The Map specialization suitable for this CrsMatrix specialization.
    typedef Map<LocalOrdinal,GlobalOrdinal,Node>  map_type;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Constructor specifying fixed number of entries for each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
               size_t maxNumEntriesPerRow,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying (possibly different) number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
               const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and fixed number of entries for each row.
    ///
    /// The column Map will be used to filter any matrix entries
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
               const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
               size_t maxNumEntriesPerRow,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// The column Map will be used to filter any matrix indices
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
               const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
               const ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a previously constructed graph.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrix with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graph must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrix constructor.
    ///
    /// This constructor is marked \c explicit so that you can't
    /// create a CrsMatrix by accident when passing a CrsGraph into a
    /// function that takes a CrsMatrix.
    ///
    /// \param graph [in] The graph structure of the sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit CrsMatrix (const Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >& graph,
                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and arrays containing the matrix in sorted, local ids.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param rowPointers [in] The beginning of each row in the matrix,
    ///   as in a CSR "rowptr" array.  The length of this vector should be
    ///   equal to the number of rows in the graph, plus one.  This last
    ///   entry should store the nunber of nonzeros in the matrix.
    ///
    /// \param columnIndices [in] The local indices of the columns,
    ///   as in a CSR "colind" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param values [in] The local entries in the matrix,
    ///   as in a CSR "vals" array.  The length of this vector
    ///   should be equal to the number of unknowns in the matrix.
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrix (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
               const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
               const ArrayRCP<size_t> & rowPointers,
               const ArrayRCP<LocalOrdinal> & columnIndices,
               const ArrayRCP<Scalar> & values,
               const RCP<ParameterList>& params = null);

    // This friend declaration makes the clone() method work.
    template <class S2, class LO2, class GO2, class N2, class LMO2>
    friend class CrsMatrix;

    /// \brief Create a deep copy of this CrsMatrix, where the copy
    ///   may have a different Node type.
    ///
    /// \param node2 [in] Kokkos Node instance for the returned copy.
    /// \param params [in/out] Optional list of parameters. If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    ///
    /// Parameters to \c params:
    /// - "Static profile clone" [boolean, default: true] If \c true,
    ///   create the copy with a static allocation profile. If false,
    ///   use a dynamic allocation profile.
    /// - "Locally indexed clone" [boolean] If \c true, fill clone
    ///   using this matrix's column Map and local indices.  This
    ///   matrix must have a column Map in order for this to work.  If
    ///   false, fill clone using global indices.  By default, this
    ///   will use local indices only if this matrix is using local
    ///   indices.
    /// - "fillComplete clone" [boolean, default: true] If \c true,
    ///   call fillComplete() on the cloned CrsMatrix object, with
    ///   parameters from the input parameters' "CrsMatrix" sublist
    ///   The domain Map and range Map passed to fillComplete() are
    ///   those of the map being cloned, if they exist. Otherwise, the
    ///   row Map is used.
    template <class Node2>
    RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node2, typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Node2>::SparseOps> >
    clone (const RCP<Node2>& node2, const RCP<ParameterList>& params = null) const
    {
      const char tfecfFuncName[] = "clone";

      // Get parameter values.  Set them initially to their default values.
      bool fillCompleteClone = true;
      bool useLocalIndices = this->hasColMap ();
      ProfileType pftype = StaticProfile;
      if (! params.is_null ()) {
        fillCompleteClone = params->get ("fillComplete clone", fillCompleteClone);
        useLocalIndices = params->get ("Locally indexed clone", useLocalIndices);

        bool staticProfileClone = true;
        staticProfileClone = params->get ("Static profile clone", staticProfileClone);
        pftype = staticProfileClone ? StaticProfile : DynamicProfile;
      }

      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        ! this->hasColMap () && useLocalIndices, std::runtime_error,
        ": You requested that the returned clone have local indices, but the "
        "the source matrix does not have a column Map yet.");

      typedef typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Node2>::SparseOps LocalMatOps2;
      typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node2, LocalMatOps2> CrsMatrix2;
      typedef Map<LocalOrdinal, GlobalOrdinal, Node2> Map2;
      RCP<const Map2> clonedRowMap = this->getRowMap ()->template clone<Node2> (node2);

      RCP<CrsMatrix2> clonedMatrix;
      ArrayRCP<const size_t> numEntries;
      size_t numEntriesForAll = 0;
      if (! staticGraph_->indicesAreAllocated ()) {
        if (! staticGraph_->numAllocPerRow_.is_null ()) {
          numEntries = staticGraph_->numAllocPerRow_;
        }
        else {
          numEntriesForAll = staticGraph_->numAllocForAllRows_;
        }
      }
      else if (! staticGraph_->numRowEntries_.is_null ()) {
        numEntries = staticGraph_->numRowEntries_;
      }
      else if (staticGraph_->nodeNumAllocated_ == 0) {
        numEntriesForAll = 0;
      }
      else {
        // We're left with the case that we have optimized storage.
        // In this case, we have to construct a list of row sizes.
        TEUCHOS_TEST_FOR_EXCEPTION(
          getProfileType() != StaticProfile, std::logic_error,
          "Internal logic error. Please report this to Tpetra team." )

        const size_t numRows = this->getNodeNumRows ();
        numEntriesForAll = 0;
        ArrayRCP<size_t> numEnt;
        if (numRows != 0) {
          numEnt = arcp<size_t> (numRows);
        }
        for (size_t i = 0; i < numRows; ++i) {
          numEnt[i] = staticGraph_->rowPtrs_[i+1] - staticGraph_->rowPtrs_[i];
        }
        numEntries = numEnt;
      }

      RCP<ParameterList> matrixparams =
        params.is_null () ? null : sublist (params,"CrsMatrix");
      if (useLocalIndices) {
        RCP<const Map2> clonedColMap =
          this->getColMap ()->template clone<Node2> (node2);
        if (numEntries.is_null ()) {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, clonedColMap,
                                              numEntriesForAll, pftype,
                                              matrixparams));
        }
        else {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, clonedColMap,
                                              numEntries, pftype,
                                              matrixparams));
        }
      }
      else {
        if (numEntries.is_null ()) {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, numEntriesForAll,
                                              pftype, matrixparams));
        }
        else {
          clonedMatrix = rcp (new CrsMatrix2 (clonedRowMap, numEntries, pftype,
                                              matrixparams));
        }
      }
      // done with these
      numEntries = null;
      numEntriesForAll = 0;

      if (useLocalIndices) {
        clonedMatrix->allocateValues (LocalIndices,
                                      CrsMatrix2::GraphNotYetAllocated);
        if (this->isLocallyIndexed ()) {
          ArrayView<const LocalOrdinal> linds;
          ArrayView<const Scalar>       vals;
          for (LocalOrdinal lrow = clonedRowMap->getMinLocalIndex ();
               lrow <= clonedRowMap->getMaxLocalIndex ();
               ++lrow) {
            this->getLocalRowView (lrow, linds, vals);
            if (linds.size ()) {
              clonedMatrix->insertLocalValues (lrow, linds, vals);
            }
          }
        }
        else { // this->isGloballyIndexed()
          Array<LocalOrdinal> linds;
          Array<Scalar> vals;
          for (LocalOrdinal lrow = clonedRowMap->getMinLocalIndex ();
               lrow <= clonedRowMap->getMaxLocalIndex ();
               ++lrow) {
            size_t theNumEntries = this->getNumEntriesInLocalRow (lrow);
            if (theNumEntries > Teuchos::as<size_t> (linds.size ())) {
              linds.resize (theNumEntries);
            }
            if (theNumEntries > Teuchos::as<size_t> (vals.size ())) {
              vals.resize (theNumEntries);
            }
            this->getLocalRowCopy (clonedRowMap->getGlobalElement (lrow),
                                   linds (), vals (), theNumEntries);
            if (theNumEntries != 0) {
              clonedMatrix->insertLocalValues (lrow, linds (0, theNumEntries),
                                               vals (0, theNumEntries));
            }
          }
        }
      }
      else { // useGlobalIndices
        clonedMatrix->allocateValues (GlobalIndices,
                                      CrsMatrix2::GraphNotYetAllocated);
        if (this->isGloballyIndexed ()) {
          ArrayView<const GlobalOrdinal> ginds;
          ArrayView<const Scalar>         vals;
          for (GlobalOrdinal grow = clonedRowMap->getMinGlobalIndex ();
               grow <= clonedRowMap->getMaxGlobalIndex ();
               ++grow) {
            this->getGlobalRowView (grow, ginds, vals);
            if (ginds.size () > 0) {
              clonedMatrix->insertGlobalValues (grow, ginds, vals);
            }
          }
        }
        else { // this->isLocallyIndexed()
          Array<GlobalOrdinal> ginds;
          Array<Scalar> vals;
          for (GlobalOrdinal grow = clonedRowMap->getMinGlobalIndex ();
               grow <= clonedRowMap->getMaxGlobalIndex ();
               ++grow) {
            size_t theNumEntries = this->getNumEntriesInGlobalRow (grow);
            if (theNumEntries > Teuchos::as<size_t> (ginds.size ())) {
              ginds.resize (theNumEntries);
            }
            if (theNumEntries > Teuchos::as<size_t> (vals.size ())) {
              vals.resize (theNumEntries);
            }
            this->getGlobalRowCopy (grow, ginds (), vals (), theNumEntries);
            if (theNumEntries != 0) {
              clonedMatrix->insertGlobalValues (grow, ginds (0, theNumEntries),
                                                vals (0, theNumEntries));
            }
          }
        }
      }

      if (fillCompleteClone) {
        RCP<ParameterList> fillparams =
          params.is_null () ? Teuchos::null : sublist (params, "fillComplete");
        try {
          RCP<const Map2> clonedRangeMap;
          RCP<const Map2> clonedDomainMap;
          if (! this->getRangeMap ().is_null () &&
              this->getRangeMap () != clonedRowMap) {
            clonedRangeMap  = this->getRangeMap ()->template clone<Node2> (node2);
          }
          else {
            clonedRangeMap = clonedRowMap;
          }
          if (! this->getDomainMap ().is_null () &&
              this->getDomainMap () != clonedRowMap) {
            clonedDomainMap = this->getDomainMap ()->template clone<Node2> (node2);
          }
          else {
            clonedDomainMap = clonedRowMap;
          }
          clonedMatrix->fillComplete (clonedDomainMap, clonedRangeMap,
                                      fillparams);
        }
        catch (std::exception &e) {
          const bool caughtExceptionOnClone = true;
          TEUCHOS_TEST_FOR_EXCEPTION(
            caughtExceptionOnClone, std::runtime_error,
            Teuchos::typeName (*this) << std::endl << "clone: " << std::endl <<
            "Caught the following exception while calling fillComplete() on a "
            "clone of type" << std::endl << Teuchos::typeName (*clonedMatrix)
            << ": " << std::endl << e.what () << std::endl);
        }
      }
      return clonedMatrix;
    }

    //! Destructor.
    virtual ~CrsMatrix ();

    //@}
    //! @name Methods for inserting, modifying, or removing entries
    //@{

    /// \brief Insert one or more entries into the matrix, using global indices.
    ///
    /// \param globalRow [in] Global index of the row into which to
    ///   insert the entries.
    /// \param cols [in] Global indices of the columns into which
    ///   to insert the entries.
    /// \param values [in] Values to insert into the above columns.
    ///
    /// For all k in 0, ..., <tt>col.size()-1</tt>, insert the value
    /// <tt>values[k]</tt> into entry <tt>(globalRow, cols[k])</tt> of
    /// the matrix.  If that entry already exists, add the new value
    /// to the old value.
    ///
    /// This is a local operation.  It does not communicate (using
    /// MPI).  If row \c globalRow is owned by the calling process,
    /// the entries will be inserted immediately.  Otherwise, if that
    /// row is <i>not</i> owned by the calling process, then the
    /// entries will be stored locally for now, and only communicated
    /// to the process that owns the row when either fillComplete() or
    /// globalAssemble() is called.  If that process already has an
    /// entry, the incoming value will be added to the old value, just
    /// as if it were inserted on the owning process.
    //
    /// If the matrix has a column Map (<tt>hasColMap() == true</tt>),
    /// and if globalRow is owned by process p, then it is forbidden
    /// to insert column indices that are not in the column Map on
    /// process p.  Tpetra will test the input column indices to
    /// ensure this is the case, but if \c globalRow is not owned by
    /// the calling process, the test will be deferred until the next
    /// call to globalAssemble() or fillComplete().
    ///
    /// \warning The behavior described in the above paragraph differs
    ///   from that of Epetra.  If the matrix has a column Map,
    ///   Epetra_CrsMatrix "filters" column indices not in the column
    ///   Map.  Many users found this confusing, so we changed it so
    ///   that nonowned column indices are forbidden.
    ///
    /// It is legal to call this method whether the matrix's column
    /// indices are globally or locally indexed.  If the matrix's
    /// column indices are locally indexed (<tt>isLocallyIndexed() ==
    /// true</tt>), then this method will convert the input global
    /// column indices to local column indices.
    ///
    /// For better performance when filling entries into a sparse
    /// matrix, consider the following tips:
    /// <ol>
    /// <li>Use local indices (e.g., insertLocalValues()) if you know
    ///   the column Map in advance.  Converting global indices to
    ///   local indices is expensive.  Of course, if you don't know
    ///   the column Map in advance, you must use global indices.</li>
    /// <li>When invoking the CrsMatrix constructor, give the best
    ///   possible upper bounds on the number of entries in each row
    ///   of the matrix.  This will avoid expensive reallocation if
    ///   your bound was not large enough.</li>
    /// <li>If your upper bound on the number of entries in each row
    ///   will always be correct, create the matrix with
    ///   StaticProfile.  This uses a faster and more compact data
    ///   structure to store the matrix.</li>
    /// <li>If you plan to reuse a matrix's graph structure, but
    ///   change its values, in repeated fillComplete() / resumeFill()
    ///   cycles, you can get the best performance by creating the
    ///   matrix with a const CrsGraph.  Do this by using the
    ///   CrsMatrix constructor that accepts an RCP of a const
    ///   CrsGraph.  If you do this, you must use the "replace" or
    ///   "sumInto" methods to change the values of the matrix; you
    ///   may not use insertGlobalValues() or
    ///   insertLocalValues().</li>
    /// </ol>
    void
    insertGlobalValues (const GlobalOrdinal globalRow,
                        const ArrayView<const GlobalOrdinal>& cols,
                        const ArrayView<const Scalar>& vals);

    /// \brief Insert one or more entries into the matrix, using local indices.
    ///
    /// \param LocalRow [in] Local index of the row into which to
    ///   insert the entries.  It must be owned by the row Map on the
    ///   calling process.
    /// \param cols [in] Local indices of the columns into which to
    ///   insert the entries.  All of the column indices must be owned
    ///   by the column Map on the calling process.
    /// \param values [in] Values to insert into the above columns.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, insert the value
    /// <tt>values[k]</tt> into entry <tt>(globalRow, cols[k])</tt> of
    /// the matrix.  If that entry already exists, add the new value
    /// to the old value.
    ///
    /// In order to call this method, the matrix must be locally
    /// indexed, and it must have a column Map.
    ///
    /// For better performance when filling entries into a sparse
    /// matrix, consider the following tips:
    /// <ol>
    /// <li>When invoking the CrsMatrix constructor, give the best
    ///   possible upper bounds on the number of entries in each row
    ///   of the matrix.  This will avoid expensive reallocation if
    ///   your bound was not large enough.</li>
    /// <li>If your upper bound on the number of entries in each row
    ///   will always be correct, create the matrix with
    ///   StaticProfile.  This uses a faster and more compact data
    ///   structure to store the matrix.</li>
    /// <li>If you plan to reuse a matrix's graph structure, but
    ///   change its values, in repeated fillComplete() / resumeFill()
    ///   cycles, you can get the best performance by creating the
    ///   matrix with a const CrsGraph.  Do this by using the
    ///   CrsMatrix constructor that accepts an RCP of a const
    ///   CrsGraph.  If you do this, you must use the "replace" or
    ///   "sumInto" methods to change the values of the matrix; you
    ///   may not use insertGlobalValues() or
    ///   insertLocalValues().</li>
    /// </ol>
    void
    insertLocalValues (const LocalOrdinal localRow,
                       const ArrayView<const LocalOrdinal> &cols,
                       const ArrayView<const Scalar> &vals);

    /// Replace one or more entries' values, using global indices.
    ///
    /// \param globalRow [in] Global index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param cols [in] Global indices of the columns in which to
    ///   replace the entries.
    /// \param vals [in] Values to use for replacing the entries.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, replace the value
    /// at entry <tt>(globalRow, cols[k])</tt> of the matrix with
    /// <tt>vals[k]</tt>.  That entry must exist in the matrix
    /// already.
    ///
    /// If <tt>(globalRow, cols[k])</tt> corresponds to an entry that
    /// is duplicated in this matrix row (likely because it was
    /// inserted more than once and fillComplete() has not been called
    /// in the interim), the behavior of this method is not defined.
    void
    replaceGlobalValues (GlobalOrdinal globalRow,
                         const ArrayView<const GlobalOrdinal>& cols,
                         const ArrayView<const Scalar>& vals);

    /// Replace one or more entries' values, using global indices.
    ///
    /// \param localRow [in] local index of the row in which to
    ///   replace the entries.  This row <i>must</i> be owned by the
    ///   calling process.
    /// \param cols [in] Local indices of the columns in which to
    ///   replace the entries.
    /// \param vals [in] Values to use for replacing the entries.
    ///
    /// For all k in 0, ..., <tt>cols.size()-1</tt>, replace the value
    /// at entry <tt>(localRow, cols[k])</tt> of the matrix with
    /// <tt>vals[k]</tt>.  That entry must exist in the matrix
    /// already.
    void
    replaceLocalValues (LocalOrdinal localRow,
                        const ArrayView<const LocalOrdinal> &cols,
                        const ArrayView<const Scalar>       &vals);

    /// \brief Sum into one or more sparse matrix entries, using global indices.
    ///
    /// This is a local operation; it does not involve communication.
    /// However, if you sum into rows not owned by the calling
    /// process, it may result in future communication in
    /// globalAssemble() (which is called by fillComplete()).
    ///
    /// If globalRow is owned by the calling process, then this method
    /// performs the sum-into operation right away.  Otherwise, if the
    /// row is <i>not</i> owned by the calling process, this method
    /// defers the sum-into operation until globalAssemble().  That
    /// method communicates data for nonowned rows to the processes
    /// that own those rows.  Then, globalAssemble() does one of the
    /// following:
    /// - It calls insertGlobalValues() for that data if the matrix
    ///   has a dynamic graph.
    /// - It calls sumIntoGlobalValues() for that data if the matrix
    ///   has a static graph.  The matrix silently ignores
    ///   (row,column) pairs that do not exist in the graph.
    ///
    /// \param globalRow [in] The global index of the row in which to
    ///   sum into the matrix entries.
    /// \param cols [in] One or more column indices.
    /// \param vals [in] One or more values corresponding to those
    ///   column indices.  <tt>vals[k]</tt> corresponds to
    ///   <tt>cols[k]</tt>.
    void
    sumIntoGlobalValues (const GlobalOrdinal globalRow,
                         const ArrayView<const GlobalOrdinal> &cols,
                         const ArrayView<const Scalar>        &vals);

    /// \brief Sum into one or more sparse matrix entries, using local indices.
    ///
    /// \param localRow [in] Local index of a row.  This row must be
    ///   owned by the calling process.
    /// \param cols [in] Local indices of the columns whose entries we
    ///   want to modify.
    /// \param vals [in] Values corresponding to the above column
    ///   indices.  <tt>vals[k]</tt> corresponds to <tt>cols[k]</tt>.
    void
    sumIntoLocalValues (const LocalOrdinal localRow,
                        const ArrayView<const LocalOrdinal>  &cols,
                        const ArrayView<const Scalar>        &vals);

    //! Set all matrix entries equal to \c alpha.
    void setAllToScalar (const Scalar &alpha);

    //! Scale the matrix's values: <tt>this := alpha*this</tt>.
    void scale (const Scalar &alpha);

    //! Sets the 1D pointer arrays of the graph.
    /**
       \pre <tt>hasColMap() == true</tt>
       \pre <tt>getGraph() != Teuchos::null</tt>
       \pre No insert/sum routines have been called

       \warning This method is intended for expert developer use only, and should never be called by user code.
    */
    void
    setAllValues (const ArrayRCP<size_t>& rowPointers,
                  const ArrayRCP<LocalOrdinal>& columnIndices,
                  const ArrayRCP<Scalar>& values);

    void getAllValues(ArrayRCP<const size_t> & rowPointers,ArrayRCP<const LocalOrdinal> & columnIndices, ArrayRCP<const Scalar> & values) const;


    //@}
    //! @name Transformational methods
    //@{

    /// \brief Communicate nonlocal contributions to other processes.
    ///
    /// Users do not normally need to call this method.  fillComplete
    /// always calls this method, unless you specifically tell
    /// fillComplete to do otherwise by setting its "No Nonlocal
    /// Changes" parameter to \c true.  Thus, it suffices to call
    /// fillComplete.
    ///
    /// Methods like insertGlobalValues and sumIntoGlobalValues let
    /// you add or modify entries in rows that are not owned by the
    /// calling process.  These entries are called "nonlocal
    /// contributions."  The methods that allow nonlocal contributions
    /// store the entries on the calling process, until globalAssemble
    /// is called.  globalAssemble sends these nonlocal contributions
    /// to the process(es) that own them, where they then become part
    /// of the matrix.
    ///
    /// This method only does global assembly if there are nonlocal
    /// entries on at least one process.  It does an all-reduce to
    /// find that out.  If not, it returns early, without doing any
    /// more communication or work.
    ///
    /// If you previously inserted into a row which is not owned by
    /// <i>any</i> process in the row Map, the behavior of this method
    /// is undefined.  It may detect the invalid row indices and throw
    /// an exception, or it may silently drop the entries inserted
    /// into invalid rows.  Behavior may vary, depending on whether
    /// Tpetra was built with debug checking enabled.
    void globalAssemble();

    /// \brief Resume operations that may change the values or
    ///   structure of the matrix.
    ///
    /// This method must be called as a collective operation.
    ///
    /// Calling fillComplete "freezes" both the values and the
    /// structure of the matrix.  If you want to modify the matrix
    /// again, you must first call resumeFill.  You then may not call
    /// resumeFill again on that matrix until you first call
    /// fillComplete.  You may make sequences of fillComplete,
    /// resumeFill calls as many times as you wish.
    ///
    /// \post <tt>isFillActive() && ! isFillComplete()</tt>
    void resumeFill (const RCP<ParameterList>& params = null);

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

      Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>

      Parameters:

      - "No Nonlocal Changes" (\c bool): Default is false.  If true,
        the caller promises that no modifications to nonowned rows
        have happened on any process since the last call to
        fillComplete.  This saves a global all-reduce to check whether
        any process did a nonlocal insert.  Nonlocal changes include
        any sumIntoGlobalValues or insertGlobalValues call with a row
        index that is not in the row Map of the calling process.
    */
    void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap,
                      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap,
                      const RCP<ParameterList> &params = null);

    /*! \brief Signal that data entry is complete.

      Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

      \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>
    */
    void fillComplete(const RCP<ParameterList> &params = null);

    /// \brief Perform a fillComplete on a matrix that already has data.
    ///
    /// The matrix must already have filled local 1-D storage
    /// (lclInds1D_ and rowPtrs_ for the graph, and values1D_ in the
    /// matrix).  If the matrix has been constructed in any other way,
    /// this method will throw an exception.  This routine is needed
    /// to support other Trilinos packages and should not be called by
    /// ordinary users.
    ///
    /// \warning This method is intended for expert developer use
    ///   only, and should never be called by user code.
    void
    expertStaticFillComplete (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                              const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                              const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer=Teuchos::null,
                              const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter=Teuchos::null,
                              const RCP<ParameterList> &params=Teuchos::null);

    /// \brief Replace the current domain Map and Import with the given objects.
    ///
    /// \param newDomainMap [in] New domain Map.  Must be nonnull.
    /// \param newImporter [in] Optional Import object.  If null, we
    ///   will compute it.
    ///
    /// \pre The matrix must be fill complete:
    ///   <tt>isFillComplete() == true</tt>.
    /// \pre If the Import is provided, its target Map must be the
    ///   same as the column Map of the matrix.
    /// \pre If the Import is provided, its source Map must be the
    ///   same as the provided new domain Map.
    void
    replaceDomainMapAndImporter (const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& newDomainMap,
                                 Teuchos::RCP<const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> >& newImporter);

    /// \brief Remove processes owning zero rows from the Maps and their communicator.
    ///
    /// \warning This method is ONLY for use by experts.  We highly
    ///   recommend using the nonmember function of the same name
    ///   defined in Tpetra_DistObject_decl.hpp.
    ///
    /// \warning We make NO promises of backwards compatibility.
    ///   This method may change or disappear at any time.
    ///
    /// \param newMap [in] This <i>must</i> be the result of calling
    ///   the removeEmptyProcesses() method on the row Map.  If it
    ///   is not, this method's behavior is undefined.  This pointer
    ///   will be null on excluded processes.
    virtual void
    removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap);

    //@}
    //! @name Methods implementing RowMatrix
    //@{

    //! The communicator over which the matrix is distributed.
    RCP<const Comm<int> > getComm() const;

    //! The Kokkos Node instance.
    RCP<Node> getNode() const;

    //! The Map that describes the row distribution in this matrix.
    RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRowMap() const;

    //! The Map that describes the column distribution in this matrix.
    RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getColMap() const;

    //! This matrix's graph, as a RowGraph.
    RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const;

    //! This matrix's graph, as a CrsGraph.
    RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > getCrsGraph() const;

    /// \brief Number of global elements in the row map of this matrix.
    ///
    /// This is <it>not</it> the number of rows in the matrix as a
    /// mathematical object.  This method returns the global sum of
    /// the number of local elements in the row map on each processor,
    /// which is the row map's getGlobalNumElements().  Since the row
    /// map is not one-to-one in general, that global sum could be
    /// different than the number of rows in the matrix.  If you want
    /// the number of rows in the matrix, ask the range map for its
    /// global number of elements, using the following code:
    /// <code>
    /// global_size_t globalNumRows = getRangeMap()->getGlobalNumElements();
    /// </code>
    /// This method retains the behavior of Epetra, which also asks
    /// the row map for the global number of rows, rather than asking
    /// the range map.
    ///
    /// \warning Undefined if isFillActive().
    ///
    global_size_t getGlobalNumRows() const;

    /// \brief The number of global columns in the matrix.
    ///
    /// This equals the number of entries in the matrix's domain Map.
    ///
    /// \warning Undefined if isFillActive().
    global_size_t getGlobalNumCols() const;

    /// \brief The number of matrix rows owned by the calling process.
    ///
    /// Note that the sum of all the return values over all processes
    /// in the row Map's communicator does not necessarily equal the
    /// global number of rows in the matrix, if the row Map is
    /// overlapping.
    size_t getNodeNumRows() const;

    /// \brief The number of columns connected to the locally owned rows of this matrix.
    ///
    /// Throws std::runtime_error if <tt>! hasColMap ()</tt>.
    size_t getNodeNumCols() const;

    //! The index base for global indices for this matrix.
    GlobalOrdinal getIndexBase() const;

    //! The global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const;

    //! The local number of entries in this matrix.
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

    /// \brief Whether the matrix is locally indexed on the calling process.
    ///
    /// The matrix is locally indexed on the calling process if and
    /// only if all of the following hold:
    /// <ol>
    /// <li> The matrix is not empty on the calling process </li>
    /// <li> The matrix has a column Map </li>
    /// </ol>
    ///
    /// The following is always true:
    /// \code
    /// (! locallyIndexed() && ! globallyIndexed()) || (locallyIndexed() || globallyIndexed());
    /// \endcode
    /// That is, a matrix may be neither locally nor globally indexed,
    /// but it can never be both.  Furthermore a matrix that is not
    /// fill complete, might have some processes that are neither
    /// locally nor globally indexed, and some processes that are
    /// globally indexed.  The processes that are neither do not have
    /// any entries.
    bool isLocallyIndexed() const;

    /// \brief Whether the matrix is globally indexed on the calling process.
    ///
    /// The matrix is globally indexed on the calling process if and
    /// only if all of the following hold:
    /// <ol>
    /// <li> The matrix is not empty on the calling process </li>
    /// <li> The matrix does not yet have a column Map </li>
    /// </ol>
    ///
    /// The following is always true:
    /// \code
    /// (! locallyIndexed() && ! globallyIndexed()) || (locallyIndexed() || globallyIndexed());
    /// \endcode
    /// That is, a matrix may be neither locally nor globally indexed,
    /// but it can never be both.  Furthermore a matrix that is not
    /// fill complete, might have some processes that are neither
    /// locally nor globally indexed, and some processes that are
    /// globally indexed.  The processes that are neither do not have
    /// any entries.
    bool isGloballyIndexed() const;

    /// \brief Whether the matrix is fill complete.
    ///
    /// A matrix is <i>fill complete</i> (or "in compute mode") when
    /// fillComplete() has been called without an intervening call to
    /// resumeFill().  A matrix must be fill complete in order to call
    /// computational kernels like sparse matrix-vector multiply and
    /// sparse triangular solve.  A matrix must be <i>not</i> fill
    /// complete ("in edit mode") in order to call methods that
    /// insert, modify, or remove entries.
    ///
    /// The following are always true:
    /// <ul>
    /// <li> <tt> isFillActive() == ! isFillComplete() </tt>
    /// <li> <tt> isFillActive() || isFillComplete() </tt>
    /// </ul>
    ///
    /// A matrix starts out (after its constructor returns) as not
    /// fill complete.  It becomes fill complete after fillComplete()
    /// returns, and becomes not fill complete again if resumeFill()
    /// is called.  Some methods like clone() and some of the
    /// "nonmember constructors" (like importAndFillComplete() and
    /// exportAndFillComplete()) may return a fill-complete matrix.
    bool isFillComplete() const;

    /// \brief Whether the matrix is not fill complete.
    ///
    /// A matrix is <i>fill complete</i> (or "in compute mode") when
    /// fillComplete() has been called without an intervening call to
    /// resumeFill().  A matrix must be fill complete in order to call
    /// computational kernels like sparse matrix-vector multiply and
    /// sparse triangular solve.  A matrix must be <i>not</i> fill
    /// complete ("in edit mode") in order to call methods that
    /// insert, modify, or remove entries.
    ///
    /// The following are always true:
    /// <ul>
    /// <li> <tt> isFillActive() == ! isFillComplete() </tt>
    /// <li> <tt> isFillActive() || isFillComplete() </tt>
    /// </ul>
    ///
    /// A matrix starts out (after its constructor returns) as not
    /// fill complete.  It becomes fill complete after fillComplete()
    /// returns, and becomes not fill complete again if resumeFill()
    /// is called.  Some methods like clone() and some of the
    /// "nonmember constructors" (like importAndFillComplete() and
    /// exportAndFillComplete()) may return a fill-complete matrix.
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


    //! Returns \c true if getLocalRowView() and getGlobalRowView() are valid for this class
    virtual bool supportsRowViews() const;

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
    void
    getGlobalRowCopy (GlobalOrdinal GlobalRow,
                      const ArrayView<GlobalOrdinal> &Indices,
                      const ArrayView<Scalar> &Values,
                      size_t &NumEntries) const;

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
    void
    getLocalRowCopy (LocalOrdinal LocalRow,
                     const ArrayView<LocalOrdinal> &Indices,
                     const ArrayView<Scalar> &Values,
                     size_t &NumEntries) const;

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \param Values    - (Out) Row values
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    void
    getGlobalRowView (GlobalOrdinal GlobalRow,
                      ArrayView<const GlobalOrdinal> &indices,
                      ArrayView<const Scalar> &values) const;

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \param Values   - (Out) Row values
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    void
    getLocalRowView (LocalOrdinal LocalRow,
                     ArrayView<const LocalOrdinal> &indices,
                     ArrayView<const Scalar> &values) const;

    //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
    /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
      the zero and non-zero diagonals owned by this node. */
    void getLocalDiagCopy (Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const;

    /// \brief Get offsets of the diagonal entries in the matrix.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    /// \warning This method must be called collectively.  We reserve
    ///   the right to do extra checking in a debug build that will
    ///   require collectives.
    ///
    /// \pre The matrix must be locally indexed (which means that it
    ///   has a column Map).
    /// \pre All diagonal entries of the matrix's graph must be
    ///   populated on this process.  Results are undefined otherwise.
    /// \post <tt>offsets.size() == getNodeNumRows()</tt>
    ///
    /// This method creates an array of offsets of the local diagonal
    /// entries in the matrix.  This array is suitable for use in the
    /// two-argument version of getLocalDiagCopy().  However, its
    /// contents are not defined in any other context.  For example,
    /// you should not rely on offsets[i] being the index of the
    /// diagonal entry in the views returned by getLocalRowView().
    /// This may be the case, but it need not be.  (For example, we
    /// may choose to optimize the lookups down to the optimized
    /// storage level, in which case the offsets will be computed with
    /// respect to the underlying storage format, rather than with
    /// respect to the views.)
    ///
    /// Calling any of the following invalidates the output array:
    /// - insertGlobalValues()
    /// - insertLocalValues()
    /// - fillComplete() (with a dynamic graph)
    /// - resumeFill() (with a dynamic graph)
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().  "Invalidates" means
    /// that you must call this method again to recompute the offsets.
    void getLocalDiagOffsets (Teuchos::ArrayRCP<size_t>& offsets) const;

    /// \brief Variant of getLocalDiagCopy() that uses precomputed offsets.
    ///
    /// \warning This method is only for expert users.
    /// \warning We make no promises about backwards compatibility
    ///   for this method.  It may disappear or change at any time.
    ///
    /// This method uses the offsets of the diagonal entries, as
    /// precomputed by getLocalDiagOffsets(), to speed up copying the
    /// diagonal of the matrix.  The offsets must be recomputed if any
    /// of the following methods are called:
    /// - insertGlobalValues()
    /// - insertLocalValues()
    /// - fillComplete() (with a dynamic graph)
    /// - resumeFill() (with a dynamic graph)
    ///
    /// If the matrix has a const ("static") graph, and if that graph
    /// is fill complete, then the offsets array remains valid through
    /// calls to fillComplete() and resumeFill().
    void
    getLocalDiagCopy (Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& diag,
                      const Teuchos::ArrayView<const size_t>& offsets) const;

    /** \brief . */
    void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);

    /** \brief . */
    void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);

    //@}
    //! @name Advanced templated methods
    //@{

    /// \brief Compute a sparse matrix-MultiVector product local to each process.
    ///
    /// This method computes the <i>local</i> part of <tt>Y := beta*Y
    /// + alpha*Op(A)*X</tt>, where <tt>Op(A)</tt> is either \f$A\f$,
    /// \f$A^T\f$ (the transpose), or \f$A^H\f$ (the conjugate
    /// transpose).  "The local part" means that this method does no
    /// communication between processes, even if this is necessary for
    /// correctness of the matrix-vector multiply.  Use the apply()
    /// method if you want to compute the mathematical sparse
    /// matrix-vector multiply.
    ///
    /// This method is mainly of use to Tpetra developers, though some
    /// users may find it helpful if they plan to reuse the result of
    /// doing an Import on the input MultiVector for several sparse
    /// matrix-vector multiplies with matrices that have the same
    /// column Map.
    ///
    /// When <tt>Op(A)</tt> is \f$A\f$ (<tt>trans ==
    /// Teuchos::NO_TRANS</tt>), then X's Map must be the same as the
    /// column Map of this matrix, and Y's Map must be the same as the
    /// row Map of this matrix.  We say in this case that X is
    /// "post-Imported," and Y is "pre-Exported."  When <tt>Op(A)</tt>
    /// is \f$A^T\f$ or \f$A^H\f$ (\c trans is <tt>Teuchos::TRANS</tt>
    /// or <tt>Teuchos::CONJ_TRANS</tt>, then X's Map must be the same
    /// as the row Map of this matrix, and Y's Map must be the same as
    /// the column Map of this matrix.
    ///
    /// Both X and Y must have constant stride, and they may not alias
    /// one another (that is, occupy overlapping space in memory).  We
    /// may not necessarily check for aliasing, and if we do, we will
    /// only do this in a debug build.  Aliasing X and Y may cause
    /// nondeterministically incorrect results.
    ///
    /// This method is templated on the type of entries in both the
    /// input MultiVector (\c DomainScalar) and the output MultiVector
    /// (\c RangeScalar).  Thus, this method works for MultiVector
    /// objects of arbitrary type.  However, this method only performs
    /// computation local to each MPI process.  Use
    /// CrsMatrixMultiplyOp to handle global communication (the Import
    /// and Export operations for the input resp. output MultiVector),
    /// if you have a matrix with entries of a different type than the
    /// input and output MultiVector objects.
    ///
    /// If <tt>beta == 0</tt>, this operation will enjoy overwrite
    /// semantics: Y will be overwritten with the result of the
    /// multiplication, even if it contains <tt>NaN</tt>
    /// (not-a-number) floating-point entries.  Otherwise, the
    /// multiply result will be accumulated into \c Y.
    template <class DomainScalar, class RangeScalar>
    void
    localMultiply (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp trans,
                   RangeScalar alpha,
                   RangeScalar beta) const;

    /// \brief Gauss-Seidel or SOR on \f$B = A X\f$.
    ///
    /// Apply a forward or backward sweep of Gauss-Seidel or
    /// Successive Over-Relaxation (SOR) to the linear system(s) \f$B
    /// = A X\f$.  For Gauss-Seidel, set the damping factor \c omega
    /// to 1.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in B.
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector B.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param omega [in] SOR damping factor.  omega = 1 results in
    ///   Gauss-Seidel.
    /// \param direction [in] Sweep direction: KokkosClassic::Forward or
    ///   KokkosClassic::Backward.  ("Symmetric" requires interprocess
    ///   communication (before each sweep), which is not part of the
    ///   local kernel.)
    template <class DomainScalar, class RangeScalar>
    void
    localGaussSeidel (const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                      MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                      const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                      const RangeScalar& dampingFactor,
                      const KokkosClassic::ESweepDirection direction) const;

    /// \brief Solves a linear system when the underlying matrix is triangular.
    ///
    /// X is required to be post-imported, i.e., described by the
    /// column map of the matrix. Y is required to be pre-exported,
    /// i.e., described by the row map of the matrix.
    ///
    /// This method is templated on the scalar type of MultiVector
    /// objects, allowing this method to be applied to MultiVector
    /// objects of arbitrary type. However, it is recommended that
    /// solve() not be called directly; instead, use the
    /// CrsMatrixSolveOp, as it will handle the Import/Export operations
    /// required to apply a matrix with non-trivial communication needs.
    ///
    /// Both X and Y are required to have constant stride. However,
    /// unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No
    /// run-time checking will be performed in a non-debug build.
    template <class DomainScalar, class RangeScalar>
    void
    localSolve (const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Teuchos::ETransp trans) const;

    //! Returns another CrsMatrix with the same entries, but represented as a different scalar type.
    template <class T>
    RCP<CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node> > convert() const;

    //@}
    //! @name Methods implementing Operator
    //@{

    /// \brief Compute a sparse matrix-MultiVector multiply.
    ///
    /// This method computes <tt>Y := beta*Y + alpha*Op(A)*X</tt>,
    /// where <tt>Op(A)</tt> is either \f$A\f$, \f$A^T\f$ (the
    /// transpose), or \f$A^H\f$ (the conjugate transpose).
    ///
    /// If <tt>beta == 0</tt>, this operation will enjoy overwrite
    /// semantics: Y's entries will be ignored, and Y will be
    /// overwritten with the result of the multiplication, even if it
    /// contains <tt>NaN</tt> (not-a-number) floating-point entries.
    void
    apply (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
           MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>&Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = ScalarTraits<Scalar>::one(),
           Scalar beta = ScalarTraits<Scalar>::zero()) const;

    //! Whether apply() allows applying the transpose or conjugate transpose.
    bool hasTransposeApply() const;

    /// \brief The domain Map of this operator.
    ///
    /// This is \c null until fillComplete() has been called.
    RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const;

    /// \brief The range Map of this operator.
    ///
    /// This is \c null until fillComplete() has been called.
    RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const;

    //@}
    //! @name Other "apply"-like methods
    //@{

    /// \brief "Hybrid" Jacobi + (Gauss-Seidel or SOR) on \f$B = A X\f$.
    ///
    /// "Hybrid" means Successive Over-Relaxation (SOR) or
    /// Gauss-Seidel within an (MPI) process, but Jacobi between
    /// processes.  Gauss-Seidel is a special case of SOR, where the
    /// damping factor is one.
    ///
    /// The Forward or Backward sweep directions have their usual SOR
    /// meaning within the process.  Interprocess communication occurs
    /// once before the sweep, as it normally would in Jacobi.
    ///
    /// The Symmetric sweep option means two sweeps: first Forward,
    /// then Backward.  Interprocess communication occurs before each
    /// sweep, as in Jacobi.  Thus, Symmetric results in two
    /// interprocess communication steps.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its Backward
    ///   sweep) as one.
    ///
    /// \section Tpetra_CrsMatrix_gaussSeidel_Details Requirements
    ///
    /// This method has the following requirements:
    ///
    /// 1. X is in the domain Map of the matrix.
    /// 2. The domain and row Maps of the matrix are the same.
    /// 3. The column Map contains the domain Map, and both start at the same place.
    /// 4. The row Map is uniquely owned.
    /// 5. D is in the row Map of the matrix.
    /// 6. X is actually a view of a column Map multivector.
    /// 7. Neither B nor D alias X.
    ///
    /// #1 is just the usual requirement for operators: the input
    /// multivector must always be in the domain Map.  The
    /// Gauss-Seidel kernel imposes additional requirements, since it
    ///
    /// - overwrites the input multivector with the output (which
    ///   implies #2), and
    /// - uses the same local indices for the input and output
    ///   multivector (which implies #2 and #3).
    ///
    /// #3 is reasonable if the matrix constructed the column Map,
    /// because the method that does this (CrsGraph::makeColMap) puts
    /// the local GIDs (those in the domain Map) in front and the
    /// remote GIDs (not in the domain Map) at the end of the column
    /// Map.  However, if you constructed the column Map yourself, you
    /// are responsible for maintaining this invariant.  #6 lets us do
    /// the Import from the domain Map to the column Map in place.
    ///
    /// The Gauss-Seidel kernel also assumes that each process has the
    /// entire value (not a partial value to sum) of all the diagonal
    /// elements in the rows in its row Map.  (We guarantee this anyway
    /// though the separate D vector.)  This is because each element of
    /// the output multivector depends nonlinearly on the diagonal
    /// elements.  Shared ownership of off-diagonal elements would
    /// produce different results.
    void
    gaussSeidel (const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                 const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                 const Scalar& dampingFactor,
                 const ESweepDirection direction,
                 const int numSweeps) const;

    /// \brief Version of gaussSeidel(), with fewer requirements on X.
    ///
    /// This method is just like gaussSeidel(), except that X need
    /// only be in the domain Map.  This method does not require that
    /// X be a domain Map view of a column Map multivector.  As a
    /// result, this method must copy X into a domain Map multivector
    /// before operating on it.
    ///
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param B [in] Right-hand side(s), in the range Map.
    /// \param D [in] Inverse of diagonal entries of the matrix,
    ///   in the row Map.
    /// \param dampingFactor [in] SOR damping factor.  A damping
    ///   factor of one results in Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward, Backward, or
    ///   Symmetric.
    /// \param numSweeps [in] Number of sweeps.  We count each
    ///   Symmetric sweep (including both its Forward and its
    ///   Backward sweep) as one.
    /// \param zeroInitialGuess [in] If true, this method will fill X
    ///   with zeros initially.  If false, this method will assume
    ///   that X contains a possibly nonzero initial guess on input.
    ///   Note that a nonzero initial guess may impose an additional
    ///   nontrivial communication cost (an additional Import).
    ///
    /// \pre Domain, range, and row Maps of the sparse matrix are all the same.
    /// \pre No other argument aliases X.
    void
    gaussSeidelCopy (MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B,
                     const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &D,
                     const Scalar& dampingFactor,
                     const ESweepDirection direction,
                     const int numSweeps,
                     const bool zeroInitialGuess) const;

    /// \brief Implementation of RowMatrix::add: return <tt>alpha*A + beta*this</tt>.
    ///
    /// This override of the default implementation ensures that, when
    /// called on a CrsMatrix, this method always returns a CrsMatrix
    /// of exactly the same type as <tt>*this</tt>.  "Exactly the same
    /// type" means that all the template parameters match, including
    /// the fifth template parameter.  The input matrix A need not
    /// necessarily be a CrsMatrix or a CrsMatrix of the same type as
    /// <tt>*this</tt>, though this method may be able to optimize
    /// further in that case.
    virtual Teuchos::RCP<RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    add (const Scalar& alpha,
         const RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
         const Scalar& beta,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& domainMap,
         const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& rangeMap,
         const Teuchos::RCP<Teuchos::ParameterList>& params) const;

    //@}
    //! @name Implementation of Teuchos::Describable interface
    //@{

    //! A simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    void
    describe (Teuchos::FancyOStream &out,
              const Teuchos::EVerbosityLevel verbLevel =
              Teuchos::Describable::verbLevel_default) const;

    //@}
    //! @name Implementation of DistObject interface
    //@{

    virtual bool
    checkSizes (const SrcDistObject& source);

    virtual void
    copyAndPermute (const SrcDistObject& source,
                    size_t numSameIDs,
                    const ArrayView<const LocalOrdinal> &permuteToLIDs,
                    const ArrayView<const LocalOrdinal> &permuteFromLIDs);

    virtual void
    packAndPrepare (const SrcDistObject& source,
                    const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
                    Teuchos::Array<char>& exports,
                    const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                    size_t& constantNumPackets,
                    Distributor& distor);

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
    //! @name Implementation of Packable interface
    //@{

    /// \brief Pack this object's data for an Import or Export.
    ///
    /// \warning To be called only by the packAndPrepare method of
    ///   appropriate classes of DistObject.
    ///
    /// Subclasses may override this method to speed up or otherwise
    /// improve the implementation by exploiting more specific details
    /// of the subclass.
    virtual void
    pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
          Teuchos::Array<char>& exports,
          const Teuchos::ArrayView<size_t>& numPacketsPerLID,
          size_t& constantNumPackets,
          Distributor& distor) const;

    //@}
  private:
    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Import<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& importer,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    // Friend declaration for nonmember function.
    template<class CrsMatrixType>
    friend Teuchos::RCP<CrsMatrixType>
    exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                    const Export<typename CrsMatrixType::local_ordinal_type,
                                                 typename CrsMatrixType::global_ordinal_type,
                                                 typename CrsMatrixType::node_type>& exporter,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& domainMap,
                                    const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                                 typename CrsMatrixType::global_ordinal_type,
                                                                 typename CrsMatrixType::node_type> >& rangeMap,
                                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    /// \brief Import from <tt>this</tt> to the result, and fillComplete the result.
    ///
    /// This method implements the nonmember "constructor"
    /// importAndFillCompleteCrsMatrix.  It's convenient to put that
    /// function's implementation inside the CrsMatrix class, so that
    /// we don't have to put much code in the _decl header file.
    Teuchos::RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >
    importAndFillComplete (const Import<LocalOrdinal, GlobalOrdinal, Node>& importer,
                           const Teuchos::RCP<const map_type>& domainMap,
                           const Teuchos::RCP<const map_type>& rangeMap,
                           const Teuchos::RCP<Teuchos::ParameterList>& params) const;

    /// \brief Export from <tt>this</tt> to the result, and fillComplete the result.
    ///
    /// This method implements the nonmember "constructor"
    /// exportAndFillCompleteCrsMatrix.  It's convenient to put that
    /// function's implementation inside the CrsMatrix class, so that
    /// we don't have to put much code in the _decl header file.
    Teuchos::RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >
    exportAndFillComplete (const Export<LocalOrdinal, GlobalOrdinal, Node>& exporter,
                           const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
                           const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
                           const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) const;

  private:
    // We forbid copy construction by declaring this method private
    // and not implementing it.
    CrsMatrix (const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &rhs);

    // We forbid assignment (operator=) by declaring this method
    // private and not implementing it.
    CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>&
    operator= (const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &rhs);

    /// \brief Like insertGlobalValues(), but with column filtering.
    ///
    /// "Column filtering" means that if the matrix has a column Map,
    /// then this method ignores entries in columns that are not in
    /// the column Map.
    void
    insertGlobalValuesFiltered (const GlobalOrdinal globalRow,
                                const ArrayView<const GlobalOrdinal> &indices,
                                const ArrayView<const Scalar>        &values);

    /// \brief Like insertLocalValues(), but with column filtering.
    ///
    /// "Column filtering" means that if the matrix has a column Map,
    /// then this method ignores entries in columns that are not in
    /// the column Map.
    void
    insertLocalValuesFiltered (const LocalOrdinal localRow,
                               const ArrayView<const LocalOrdinal> &indices,
                               const ArrayView<const Scalar>       &values);

    /// \brief Combine in the data using the given combine mode.
    ///
    /// The copyAndPermute() and unpackAndCombine() methods use this
    /// function to combine incoming entries from the source matrix
    /// with the target matrix's current data.  This method's behavior
    /// depends on whether the target matrix (that is, this matrix)
    /// has a static graph.
    void
    combineGlobalValues (const GlobalOrdinal globalRowIndex,
                         const Teuchos::ArrayView<const GlobalOrdinal> columnIndices,
                         const Teuchos::ArrayView<const Scalar> values,
                         const Tpetra::CombineMode combineMode);

    /// \brief Transform CrsMatrix entries, using local indices.
    ///
    /// For every entry \f$A(i,j)\f$ to transform, if \f$v_{ij}\f$ is
    /// the corresponding entry of the \c values array, then we apply
    /// the binary function f to \f$A(i,j)\f$ as follows:
    /// \f[
    ///   A(i,j) := f(A(i,j), v_{ij}).
    /// \f]
    /// For example, BinaryFunction = std::plus<Scalar> does the same
    /// thing as sumIntoLocalValues(), and BinaryFunction =
    /// project2nd<Scalar,Scalar> does the same thing as
    /// replaceLocalValues().
    ///
    /// \tparam BinaryFunction The type of binary function to apply.
    ///   std::binary_function is a model for this.
    ///
    /// \pre The matrix must have a column Map.
    ///
    /// \param localRow [in] (Local) index of the row to modify.
    ///   This row <i>must</i> be owned by the calling process.
    ///
    /// \param indices [in] (Local) indices in the row to modify.
    ///   Indices not in the column Map and their corresponding values
    ///   will be ignored.
    ///
    /// \param values [in] Values to use for modification.
    ///
    /// This method works whether indices are local or global.
    /// However, it will cost more if indices are global, since it
    /// will have to convert the local indices to global indices in
    /// that case.
    template<class BinaryFunction>
    void
    transformLocalValues (LocalOrdinal localRow,
                          const Teuchos::ArrayView<const LocalOrdinal>& indices,
                          const Teuchos::ArrayView<const Scalar>        & values,
                          BinaryFunction f)
    {
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;
      typedef Node NT;
      using Teuchos::Array;
      using Teuchos::ArrayView;

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isFillActive (),
        std::runtime_error,
        "Tpetra::CrsMatrix::transformLocalValues: Fill must be active in order "
        "to call this method.  That is, isFillActive() must return true.  If "
        "you have already called fillComplete(), you need to call resumeFill() "
        "before you can replace values.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        values.size () != indices.size (),
        std::runtime_error,
        "Tpetra::CrsMatrix::transformLocalValues: values.size () = "
        << values.size () << " != indices.size () = " << indices.size ()
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! this->hasColMap (),
        std::runtime_error,
        "Tpetra::CrsMatrix::transformLocalValues: We cannot transform local "
        "indices without a column map.");
      const bool isLocalRow = getRowMap ()->isNodeLocalElement (localRow);
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isLocalRow,
        std::runtime_error,
        "Tpetra::CrsMatrix::transformLocalValues: The specified local row "
        << localRow << " does not belong to this process "
        << getRowMap ()->getComm ()->getRank () << ".");

      RowInfo rowInfo = staticGraph_->getRowInfo (localRow);
      if (indices.size () > 0) {
        ArrayView<Scalar> curVals = this->getViewNonConst (rowInfo);
        if (isLocallyIndexed ()) {
          staticGraph_->template transformLocalValues<Scalar, BinaryFunction> (rowInfo, curVals, indices, values, f);
        }
        else if (isGloballyIndexed ()) {
          // Convert the given local indices to global indices.
          const Map<LO, GO, NT>& colMap = * (this->getColMap ());
          Array<GO> gindices (indices.size ());
          typename ArrayView<const LO>::iterator lindit = indices.begin();
          typename Array<GO>::iterator           gindit = gindices.begin();
          while (lindit != indices.end()) {
            // There is no need to filter out indices not in the column
            // Map.  Those that aren't will be mapped to invalid(),
            // which transformGlobalValues() will ignore.
            *gindit++ = colMap.getGlobalElement (*lindit++);
          }
          staticGraph_->template transformGlobalValues<Scalar, BinaryFunction> (rowInfo, curVals, gindices (), values, f);
        }
      }
    }


    /// \brief Transform CrsMatrix entries, using global indices.
    ///
    /// For every entry \f$A(i,j)\f$ to transform, if \f$v_{ij}\f$ is
    /// the corresponding entry of the \c values array, then we apply
    /// the binary function f to \f$A(i,j)\f$ as follows:
    /// \f[
    ///   A(i,j) := f(A(i,j), v_{ij}).
    /// \f]
    /// For example, BinaryFunction = std::plus<Scalar> does the same
    /// thing as sumIntoLocalValues(), and BinaryFunction =
    /// project2nd<Scalar,Scalar> does the same thing as
    /// replaceLocalValues().
    ///
    /// \tparam BinaryFunction The type of binary function to apply.
    ///   std::binary_function is a model for this.
    ///
    /// \param globalRow [in] (Global) index of the row to modify.
    ///   This row <i>must</t> be owned by the calling process.
    ///
    /// \param indices [in] (Global) indices in the row to modify.
    ///   Indices not in the column Map (if the matrix already has a
    ///   column Map) and their corresponding values will be ignored.
    ///
    /// \param values [in] Values to use for modification.
    ///
    /// This method works whether indices are local or global.
    /// However, it will cost more if indices are local, since it will
    /// have to convert the local indices to global indices in that
    /// case.
    template<class BinaryFunction>
    void
    transformGlobalValues (GlobalOrdinal globalRow,
                           const Teuchos::ArrayView<const GlobalOrdinal>& indices,
                           const Teuchos::ArrayView<const Scalar>        & values,
                           BinaryFunction f)
    {
      typedef LocalOrdinal LO;
      typedef GlobalOrdinal GO;
      typedef Node NT;
      using Teuchos::Array;
      using Teuchos::ArrayView;

      TEUCHOS_TEST_FOR_EXCEPTION(
        ! isFillActive (),
        std::runtime_error,
        "Tpetra::CrsMatrix::transformGlobalValues: Fill must be active in order "
        "to call this method.  That is, isFillActive() must return true.  If "
        "you have already called fillComplete(), you need to call resumeFill() "
        "before you can replace values.");
      TEUCHOS_TEST_FOR_EXCEPTION(
        values.size () != indices.size (),
        std::runtime_error,
        "Tpetra::CrsMatrix::transformGlobalValues: values.size () = "
        << values.size () << " != indices.size () = " << indices.size ()
        << ".");

      const LO lrow = this->getRowMap()->getLocalElement(globalRow);

      if (lrow == OTL::invalid()) {
        // FIXME (mfh 16 May 2013) We're using this exception to do
        // sumIntoGlobalValues for nonowned rows, so we might want to
        // avoid the overhead of constructing the fancy exception
        // message each time if we don't plan to use it.

        // The exception test macro doesn't let you pass an additional
        // argument to the exception's constructor, so we don't use it.
        std::ostringstream os;
        os << "transformGlobalValues: The given global row index "
           << globalRow << " is not owned by the calling process (rank "
           << this->getRowMap()->getComm()->getRank() << ").";
        throw Details::InvalidGlobalRowIndex<GO> (os.str (), globalRow);
      }

      RowInfo rowInfo = staticGraph_->getRowInfo (lrow);
      if (indices.size () > 0) {
        ArrayView<Scalar> curVals = this->getViewNonConst (rowInfo);
        if (isLocallyIndexed ()) {
          // Convert global indices to local indices.
          const Map<LO, GO, NT> &colMap = * (this->getColMap ());
          Array<LO> lindices (indices.size ());
          typename ArrayView<const GO>::iterator gindit = indices.begin();
          typename Array<LO>::iterator           lindit = lindices.begin();
          while (gindit != indices.end()) {
            // There is no need to filter out indices not in the column
            // Map.  Those that aren't will be mapped to invalid(),
            // which transformLocalValues() will ignore.
            *lindit++ = colMap.getLocalElement (*gindit++);
          }
          staticGraph_->template transformLocalValues<Scalar, BinaryFunction> (rowInfo, curVals, lindices (), values, f);
        }
        else if (isGloballyIndexed ()) {
          staticGraph_->template transformGlobalValues<Scalar, BinaryFunction> (rowInfo, curVals, indices, values, f);
        }
      }
    }


  protected:
    // useful typedefs
    typedef OrdinalTraits<LocalOrdinal>                     OTL;
    typedef ScalarTraits<Scalar>                            STS;
    typedef typename STS::magnitudeType               Magnitude;
    typedef ScalarTraits<Magnitude>                         STM;
    typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>       MV;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>             V;
    typedef CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>  Graph;
    typedef typename LocalMatOps::template bind_scalar<Scalar>::other_type                    sparse_ops_type;
    typedef typename sparse_ops_type::template graph<LocalOrdinal,Node>::graph_type          local_graph_type;
    typedef typename sparse_ops_type::template matrix<Scalar,LocalOrdinal,Node>::matrix_type local_matrix_type;
    typedef Kokkos::CrsMatrix<Scalar,LocalOrdinal,typename Node::device_type> k_local_matrix_type;

    typedef Export<LocalOrdinal, GlobalOrdinal, Node> export_type;
    typedef Import<LocalOrdinal, GlobalOrdinal, Node> import_type;

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

    /// \brief Sort the entries of each row by their column indices.
    ///
    /// This only does anything if the graph isn't already sorted
    /// (i.e., ! myGraph_->isSorted ()).  This method is called in
    /// fillComplete().
    void sortEntries();

    /// \brief Merge entries in each row with the same column indices.
    ///
    /// This only does anything if the graph isn't already merged
    /// (i.e., ! myGraph_->isMerged ()).  This method is called in
    /// fillComplete().
    void mergeRedundantEntries();

    /// \brief Clear matrix properties that require collectives.
    ///
    /// This clears whatever computeGlobalConstants() (which see)
    /// computed, in preparation for changes to the matrix.  The
    /// current implementation of this method does nothing.
    ///
    /// This method is called in resumeFill().
    void clearGlobalConstants();

    /// \brief Compute matrix properties that require collectives.
    ///
    /// The corresponding Epetra_CrsGraph method computes things
    /// like the global number of nonzero entries, that require
    /// collectives over the matrix's communicator.  The current
    /// Tpetra implementation of this method does nothing.
    ///
    /// This method is called in fillComplete().
    void computeGlobalConstants();

    /// \brief Column Map MultiVector used in apply() and gaussSeidel().
    ///
    /// This is a column Map MultiVector.  It is used as the target of
    /// the forward mode Import operation (if necessary) in apply()
    /// and gaussSeidel(), and the source of the reverse mode Export
    /// operation (if necessary) in these methods.  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the methods in
    /// question are const, yet want to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MV> importMV_;

    /// \brief Row Map MultiVector used in apply().
    ///
    /// This is a row Map MultiVector.  It is uses as the source of
    /// the forward mode Export operation (if necessary) in apply()
    /// and gaussSeidel(), and the target of the reverse mode Import
    /// operation (if necessary) in these methods.  Both of these
    /// methods create this MultiVector on demand if needed, and reuse
    /// it (if possible) for subsequent calls.
    ///
    /// This is declared <tt>mutable</tt> because the methods in
    /// question are const, yet want to cache the MultiVector for
    /// later use.
    mutable Teuchos::RCP<MV> exportMV_;

    /// \brief Create a (or fetch a cached) column Map MultiVector.
    ///
    /// \param X_domainMap [in] A domain Map Multivector.  The
    ///   returned MultiVector, if nonnull, will have the same number
    ///   of columns as Y_domainMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the domain Map and the
    /// column Map are the same (so that normally we wouldn't need the
    /// column Map MultiVector), but the following (for example)
    /// holds:
    ///
    /// 1. The kernel needs a constant stride input MultiVector, but
    ///    the given input MultiVector is not constant stride.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getColumnMapMultiVector (const MV& X_domainMap,
                             const bool force = false) const;

    /// \brief Create a (or fetch a cached) row Map MultiVector.
    ///
    /// \param Y_rangeMap [in] A range Map Multivector.  The returned
    ///   MultiVector, if nonnull, will have the same number of
    ///   columns as Y_rangeMap.
    ///
    /// \param force [in] Force creating the MultiVector if it hasn't
    ///   been created already.
    ///
    /// The \c force parameter is helpful when the range Map and the
    /// row Map are the same (so that normally we wouldn't need the
    /// row Map MultiVector), but one of the following holds:
    ///
    /// 1. The kernel needs a constant stride output MultiVector,
    ///    but the given output MultiVector is not constant stride.
    ///
    /// 2. The kernel does not permit aliasing of its input and output
    ///    MultiVector arguments, but they do alias each other.
    ///
    /// We don't test for the above in this method, because it depends
    /// on the specific kernel.
    Teuchos::RCP<MV>
    getRowMapMultiVector (const MV& Y_rangeMap,
                          const bool force = false) const;

    //! Special case of apply() for <tt>mode == Teuchos::NO_TRANS</tt>.
    void
    applyNonTranspose (const MV& X_in,
                       MV& Y_in,
                       Scalar alpha,
                       Scalar beta) const;

    //! Special case of apply() for <tt>mode != Teuchos::NO_TRANS</tt>.
    void
    applyTranspose (const MV& X_in,
                    MV& Y_in,
                    const Teuchos::ETransp mode,
                    Scalar alpha,
                    Scalar beta) const;

    // matrix data accessors
    ArrayView<const Scalar>    getView(RowInfo rowinfo) const;
    ArrayView<      Scalar>    getViewNonConst(RowInfo rowinfo);
    // local Kokkos objects

    /// \brief Fill data into the local matrix.
    ///
    /// This method is only called in fillComplete(), and it is only
    /// called if the graph's structure is already fixed (that is, if
    /// the matrix does not own the graph).
    void fillLocalMatrix(const RCP<ParameterList> &params);
    /// \brief Fill data into the local graph and matrix.
    ///
    /// This method is only called in fillComplete(), and it is only
    /// called if the graph's structure is <i>not</i> already fixed
    /// (that is, if the matrix <i>does</i> own the graph).
    void fillLocalGraphAndMatrix(const RCP<ParameterList> &params);
    // debugging
    void checkInternalState() const;

    /// \name (Global) graph pointers
    ///
    /// We keep two graph pointers in order to maintain const
    /// correctness.  myGraph_ is a graph which we create internally.
    /// Operations that change the sparsity structure also modify
    /// myGraph_.  If myGraph_ != null, then staticGraph_ == myGraph_
    /// pointerwise (we set the pointers equal to each other when we
    /// create myGraph_).  myGraph_ is only null if this CrsMatrix was
    /// created using the constructor with a const CrsGraph input
    /// argument.  In this case, staticGraph_ is set to the input
    /// CrsGraph.
    //@{
    RCP<const Graph> staticGraph_;
    RCP<      Graph>     myGraph_;
    //@}

    /// The local sparse matrix kernels, after kernel optimizations.
    ///
    /// resumeFill() sets this to null.  fillComplete() initializes
    /// this object using the local graph and matrix.
    RCP<sparse_ops_type>   lclMatOps_;
    /// The local sparse matrix, before kernel optimizations.
    ///
    /// resumeFill() sets this to null.  fillLocalGraphAndMatrix() and
    /// fillLocalMatrix() initialize this.  Once fillComplete() has
    /// initialized the sparse kernels object (lclMatOps_ above), this
    /// object is set to null again.
    RCP<local_matrix_type> lclMatrix_;

    k_local_matrix_type k_lclMatrix_;
    /// \name Sparse matrix values.
    ///
    /// values1D_ represents the values assuming "1-D" compressed
    /// sparse row storage.  values2D_ represents the values as an
    /// array of arrays, one (inner) array per row of the sparse
    /// matrix.
    ///
    /// Before allocation, both arrays are null.  After allocation,
    /// one is null.  If static allocation, then values2D_ is null.
    /// If dynamic allocation, then values1D_ is null.  The allocation
    /// always matches that of graph_, as the graph does the
    /// allocation for the matrix.
    //@{
    ArrayRCP<Scalar> values1D_;
    typedef Kokkos::View<Scalar*, typename Node::device_type> t_ValuesType;
    t_ValuesType k_values1D_;
    ArrayRCP<Array<Scalar> > values2D_;
    //@}

    //! Whether the matrix is fill complete.
    bool fillComplete_;

    /// \brief Nonlocal data added using insertGlobalValues().
    ///
    /// These data are cleared by globalAssemble(), once it finishes
    /// redistributing them to their owning processes.
    ///
    /// \note For Epetra developers: Tpetra::CrsMatrix corresponds
    ///   more to Epetra_FECrsMatrix than to Epetra_CrsMatrix.  The
    ///   insertGlobalValues() method in Tpetra::CrsMatrix, unlike
    ///   its corresponding method in Epetra_CrsMatrix, allows
    ///   insertion into rows which are not owned by the calling
    ///   process.  The globalAssemble() method redistributes these
    ///   to their owning processes.
    std::map<GlobalOrdinal, Array<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

    /// \brief Cached Frobenius norm of the (global) matrix.
    ///
    /// The value -1 (in general,
    /// <tt>-Teuchos::ScalarTraits<Magnitude>::one()</tt>) means that
    /// the norm has not yet been computed, or that the values in the
    /// matrix may have changed and the norm must be recomputed.
    mutable Magnitude frobNorm_;
public:
    template<class ViewType, class OffsetViewType >
    struct pack_functor {
      typedef typename ViewType::device_type device_type;
      ViewType src;
      ViewType dest;
      OffsetViewType src_offset;
      OffsetViewType dest_offset;
      typedef typename OffsetViewType::non_const_scalar_type ScalarIndx;

      pack_functor(ViewType dest_, ViewType src_, OffsetViewType dest_offset_, OffsetViewType src_offset_):
        src(src_),dest(dest_),src_offset(src_offset_),dest_offset(dest_offset_) {};

      KOKKOS_INLINE_FUNCTION
      void operator() (size_t row) const {
        ScalarIndx i = src_offset(row);
        ScalarIndx j = dest_offset(row);
        const ScalarIndx k = dest_offset(row+1);
        for(;j<k;j++,i++) {
          dest(j) = src(i);
        }
      }
    };
  }; // class CrsMatrix

  /** \brief Non-member function to create an empty CrsMatrix given a row map and a non-zero profile.

      \return A dynamically allocated (DynamicProfile) matrix with specified number of nonzeros per row (defaults to zero).

      \relatesalso CrsMatrix
   */
 /* template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createCrsMatrix (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
                   size_t maxNumEntriesPerRow = 0,
                   const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    using Teuchos::rcp;
    typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
    return rcp (new matrix_type (map, maxNumEntriesPerRow, DynamicProfile, params));
  }

  /// \brief Nonmember CrsMatrix constructor that fuses Import and fillComplete().
  /// \relatesalso CrsMatrix
  /// \tparam CrsMatrixType A specialization of CrsMatrix.
  ///
  /// A common use case is to create an empty destination CrsMatrix,
  /// redistribute from a source CrsMatrix (by an Import or Export
  /// operation), then call fillComplete() on the destination
  /// CrsMatrix.  This constructor fuses these three cases, for an
  /// Import redistribution.
  ///
  /// Fusing redistribution and fillComplete() exposes potential
  /// optimizations.  For example, it may make constructing the column
  /// Map faster, and it may avoid intermediate unoptimized storage in
  /// the destination CrsMatrix.  These optimizations may improve
  /// performance for specialized kernels like sparse matrix-matrix
  /// multiply, as well as for redistributing data after doing load
  /// balancing.
  ///
  /// The resulting matrix is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source matrix, and its range Map is the range Map of
  /// the source matrix.
  ///
  /// \warning If the target Map of the Import is a subset of the
  ///   source Map of the Import, then you cannot use the default
  ///   range Map.  You should instead construct a nonoverlapping
  ///   version of the target Map and supply that as the nondefault
  ///   value of the range Map.
  ///
  /// \param sourceMatrix [in] The source matrix from which to
  ///   import.  The source of an Import must have a nonoverlapping
  ///   distribution.
  ///
  /// \param importer [in] The Import instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Import must be the same as the rowMap of sourceMatrix unless
  ///   the "Reverse Mode" option on the params list, in which case
  ///   the targetMap of Import must match the rowMap of the sourceMatrix
  ///
  /// \param domainMap [in] Domain Map of the returned matrix.  If
  ///   null, we use the default, which is the domain Map of the
  ///   source matrix.
  ///
  /// \param rangeMap [in] Range Map of the returned matrix.  If
  ///   null, we use the default, which is the range Map of the
  ///   source matrix.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  importAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Import<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& importer,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap = Teuchos::null,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap = Teuchos::null,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    return sourceMatrix->importAndFillComplete (importer, domainMap, rangeMap, params);
  }

  /// \brief Nonmember CrsMatrix constructor that fuses Export and fillComplete().
  /// \relatesalso CrsMatrix
  /// \tparam CrsMatrixType A specialization of CrsMatrix.
  ///
  /// For justification, see the documentation of
  /// importAndFillCompleteCrsMatrix() (which is the Import analog of
  /// this function).
  ///
  /// The resulting matrix is fill complete (in the sense of
  /// isFillComplete()) and has optimized storage (in the sense of
  /// isStorageOptimized()).  By default, its domain Map is the domain
  /// Map of the source matrix, and its range Map is the range Map of
  /// the source matrix.
  ///
  /// \param sourceMatrix [in] The source matrix from which to
  ///   export.  Its row Map may be overlapping, since the source of
  ///   an Export may be overlapping.
  ///
  /// \param exporter [in] The Export instance containing a
  ///   precomputed redistribution plan.  The source Map of the
  ///   Export must be the same as the row Map of sourceMatrix.
  ///
  /// \param domainMap [in] Domain Map of the returned matrix.  If
  ///   null, we use the default, which is the domain Map of the
  ///   source matrix.
  ///
  /// \param rangeMap [in] Range Map of the returned matrix.  If
  ///   null, we use the default, which is the range Map of the
  ///   source matrix.
  ///
  /// \param params [in/out] Optional list of parameters.  If not
  ///   null, any missing parameters will be filled in with their
  ///   default values.
  template<class CrsMatrixType>
  Teuchos::RCP<CrsMatrixType>
  exportAndFillCompleteCrsMatrix (const Teuchos::RCP<const CrsMatrixType>& sourceMatrix,
                                  const Export<typename CrsMatrixType::local_ordinal_type,
                                               typename CrsMatrixType::global_ordinal_type,
                                               typename CrsMatrixType::node_type>& exporter,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& domainMap = Teuchos::null,
                                  const Teuchos::RCP<const Map<typename CrsMatrixType::local_ordinal_type,
                                                               typename CrsMatrixType::global_ordinal_type,
                                                               typename CrsMatrixType::node_type> >& rangeMap = Teuchos::null,
                                  const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
  {
    return sourceMatrix->exportAndFillComplete (exporter, domainMap, rangeMap, params);
  }*/
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
