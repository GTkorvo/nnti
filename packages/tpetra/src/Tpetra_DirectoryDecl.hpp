// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_DIRECTORY_DECL_HPP
#define TPETRA_DIRECTORY_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>
#include "Tpetra_MapDecl.hpp"

namespace Tpetra {

  //! Tpetra::Directory
  
  /*! For Map objects, a Directory object must be created to allow referencing
      of non-local elements. Tpetra::Directory produces and contains a uniform linear
      Map, allowing the directory to be distributed across all nodes.
      
      This class has a single constructor, accepting the Map object for which the directory 
      is created.
  */

  template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal>
  class Directory : public Teuchos::Describable {
  public:
    
    //! @name Constructors/Destructor.
    //@{ 
    
    //! Constructor
    Directory(const Map<LocalOrdinal,GlobalOrdinal> & map);
    
 public:
    //! Destructor.
    ~Directory();
    
    //@}
    
    //! @name Query methods.
    //@{ 
    
    //! \brief Returns node info for non-local Map entries.
    /*! Given a list of global IDs, this function returns the corresponding list of
      owning node IDs.

      \param globalIDs [in] List of global IDs to look up.

      \param nodeIDs [out] On return, contains node IDs for the global IDs in question. 
      -1 corresponds to global entries not present in the directory.

      \returns \c true signifies at least one specified global entry was not present in the directory; 
               \c false signifies that all specified global entries were present in the directory.

      \note If <tt>nodeIDs.size() != globalIDs.size()</tt>, then a \c std::runtime_error exception is thrown.
    */
    bool getDirectoryEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                             const Teuchos::ArrayView<int> &nodeIDs) const;
    
    //! \brief Returns node info for non-local Map entries.
    /*! Given a list of global IDs, this function returns the corresponding list of
      owning node IDs and local IDs.

      \param globalIDs [in] List of global IDs to look up.

      \param nodeIDs [out] On return, contains node IDs for the global IDs in question. 
      -1 corresponds to global entries not present in the directory.

      \param localIDs [out] On return contains the local ID of the global on the owning node. 
      Teuchos::OrdinalTraits<LocalOrdinal>::invalid() corresponds to global entries not present in the directory.

      \returns \c true signifies at least one specified global entry was not present in the directory.
               \c false signifies that all specified global entries were present in the directory.

      \note If <tt>nodeIDs.size() != globalIDs.size()</tt> or 
               <tt>localIDs.size() != globalIDs.size()</tt>, then a \c std::runtime_error exception is thrown.
    */
    bool getDirectoryEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                             const Teuchos::ArrayView<int> &nodeIDs, 
                             const Teuchos::ArrayView<LocalOrdinal> &localIDs) const;
    //@}
    
  private:
    const Map<LocalOrdinal,GlobalOrdinal> map_;
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    std::vector<GlobalOrdinal> allMinGIDs_; // size comm_->getSize()+1; entry i contains minGID for ith node, except last entry contains maxGID in the directory
    std::vector<int> nodeIDs_;
    std::vector<LocalOrdinal> LIDs_;
    Teuchos::RCP< Map<LocalOrdinal,GlobalOrdinal> > directoryMap_;

    Directory(const Directory<LocalOrdinal,GlobalOrdinal> &directory);
    
    //! declared but not defined, to prevent default implementation, do not use
    Directory<LocalOrdinal,GlobalOrdinal> & operator = (const Directory<LocalOrdinal,GlobalOrdinal> &source);

    // common code for both versions of getDirectoryEntries
    bool getEntries(const Teuchos::ArrayView<const GlobalOrdinal> &globalIDs, 
                    const Teuchos::ArrayView<int> &nodeIDs, 
                    const Teuchos::ArrayView<LocalOrdinal> &localIDs, 
                          bool computeLIDs) const;
    
    // directory setup for non-contiguous ES
    void generateDirectory();
    
  }; // class Directory
  
} // namespace Tpetra

#endif // TPETRA_DIRECTORY_DECL_HPP

