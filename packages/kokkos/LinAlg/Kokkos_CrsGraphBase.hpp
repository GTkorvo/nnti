//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#ifndef KOKKOS_CRSGRAPHBASE_HPP
#define KOKKOS_CRSGRAPHBASE_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ParameterList.hpp>
#include <numeric>

namespace Kokkos {

  /*! @class CrsGraphBase
      @brief An abstract base class providing a template for Kokkos-level sparse graph objects.
      @ingroup kokkos_crs_ops
     
      The Tpetra classes do not utilize this base class for interacting with Kokkos-level objects, and 
      the local sparse graph objects are there not required to inherit from this interface. However, 
      this class illustrates the methods that are required by Tpetra objects, and therefore provides a 
      potential starting point for extending Tpetra via new local graph and matrix types.
      
      \tparam Ordinal Defines the type of the column indices;
                      same as the LocalOrdinal template parameter of the encapsulating Tpetra::CrsGraph object.
      \tparam Node    Kokkos Node type; same as the Node template parameter of the encapsulating Tpetra::CrsGraph object.
  */
  template <class Ordinal,
            class Node>
  class CrsGraphBase {
  public:
    typedef Ordinal  ordinal_type;
    typedef Node     node_type;

    //! @name Constructors/Destructor
    //@{

    //! Default constuctor.
    CrsGraphBase(size_t numRows, const RCP<Node> &node);

    //! Destructor.
    virtual ~CrsGraphBase();

    //@}
    //! @name Accessor routines.
    //@{

    //! Node accessor.
    RCP<Node> getNode() const;

    //@}
    //! @name Data entry and accessor methods.
    //@{

    //! Return the number of rows in the graph.
    size_t getNumRows() const;

    //! Return the number of entries in the graph.
    virtual size_t getNumEntries() const = 0;

    //! Indicates that the graph has no rows or no entries.
    /**
      \note This is different from not having been finalized.
      \post If isFinalized() == false, then isEmpty() is not valid, but will \c false.
     */
    virtual bool isEmpty() const = 0;

    //! Indicatest that the graph has been finalized.
    virtual bool isFinalized() const = 0;

    //! \brief Allocate and initialize the storage for a sparse graph.
    /**
        \pre numEntries.size() == getNumRows()
     */
    virtual void allocStorage(const ArrayView<const size_t> &numEntries,
                              ArrayRCP<size_t> &ptrs, ArrayRCP<Ordinal> &inds) const;

    /// \brief Submit the indices and offset for the graph.
    /**
          \pre indices for row \c r are inds[r], where \f$r \in [b,e)\f$, where \f$b = ptrs[r]\f$ and \f$e = ptrs[r-1]\f$
          \pre ptrs has getNumRows()+1 entries
          \pre ptrs[0] == 0
          \pre ptrs[getNumRows()] == inds.size()
     */
    virtual void setStructure(const ArrayRCP<const size_t>  &ptrs,
                              const ArrayRCP<const Ordinal> &inds) = 0;

    /// \brief Finalize storage for the graph.
    ///
    /// Instruct the graph to perform any necessary manipulation.
    virtual void finalize(Teuchos::ParameterList &params);

    /* Finalize storage for the graph and matrix.
       This one cannot be made virtual without introducing a dynamic cast down to the concrete class. 
       We'll rely strictly on generic programming here.
       But it should look like this:
         void finalizeGraphAndMatrix(Matrix &mat, Teuchos::ParameterList &params);
         void finalizeMatrix(Matrix &mat, Teuchos::ParameterList &params) const;
    */

    //@}

  private:
    RCP<Node> node_;
    size_t numRows_;
  };


  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraphBase<Ordinal,Node>::CrsGraphBase(size_t numRows, const RCP<Node> &node) 
  : node_(node)
  , numRows_(numRows)
  {}

  //==============================================================================
  template <class Ordinal, class Node>
  CrsGraphBase<Ordinal,Node>::~CrsGraphBase() {
  }

  // ======= node ===========
  template <class Ordinal, class Node>
  RCP<Node> CrsGraphBase<Ordinal,Node>::getNode() const {
    return node_;
  }

  // ======= numrows ===========
  template <class Ordinal, class Node>
  size_t CrsGraphBase<Ordinal,Node>::getNumRows() const {
    return numRows_;
  }

  // ======= numrows ===========
  template <class Ordinal, class Node>
  void CrsGraphBase<Ordinal,Node>::allocStorage(
      const ArrayView<const size_t> &numEntriesPerRow,
      ArrayRCP<size_t> &ptrs, ArrayRCP<Ordinal> &inds) const
  {
    std::string tfecfFuncName("allocStorage( in(numEntriesPerRow), in(ptrs), in(inds) )");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)numEntriesPerRow.size() != (size_t)getNumRows(),
        std::runtime_error, " number of rows doesn't match.")
    const size_t totalNumEntries = std::accumulate( numEntriesPerRow.begin(), numEntriesPerRow.end(), (size_t)0 );
    // alloc inds
    if (totalNumEntries > 0) inds = arcp<Ordinal>(totalNumEntries);
    else                     inds = null;
    std::fill( inds.begin(), inds.end(), Teuchos::OrdinalTraits<Ordinal>::zero() );
    // alloc ptrs and set them
    ptrs = arcp<size_t>( getNumRows()+1 );
    ptrs[0] = 0;
    std::partial_sum( numEntriesPerRow.begin(), numEntriesPerRow.end(), ptrs.begin()+1 );
  }

  // ======= default implementation is empty ===========
  template <class Ordinal, class Node>
  void CrsGraphBase<Ordinal,Node>::finalize(Teuchos::ParameterList &params)
  { }

} // namespace Kokkos

#endif /* KOKKOS_CRSGRAPHBASE_HPP */
