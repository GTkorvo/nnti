// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef FEAPP_BLOCKDISCRETIZATION_HPP
#define FEAPP_BLOCKDISCRETIZATION_HPP

#include <vector>

#include "EpetraExt_MultiMpiComm.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockUtility.h"

#include "FEApp_AbstractDiscretization.hpp"

namespace FEApp {

  class BlockDiscretization : public FEApp::AbstractDiscretization {
  public:

    //! Constructor
    BlockDiscretization(
		 const Teuchos::RCP<const EpetraExt::MultiMpiComm>& globalComm_,
		 const Teuchos::RCP<const FEApp::AbstractDiscretization>& underlyingDisc_);

    //! Destructor
    virtual ~BlockDiscretization();

    //! Create element mesh
    virtual void createMesh();

    //! Create DOF maps
    virtual void createMaps();

    //! Create Jacobian graph
    virtual void createJacobianGraphs();

    //! Get element mesh
    virtual Teuchos::RCP<const FEApp::Mesh> 
    getMesh() const; 

    //! Get DOF map
    virtual Teuchos::RCP<const Epetra_Map> 
    getMap() const;

    //! Get overlapped DOF map
    virtual Teuchos::RCP<const Epetra_Map> 
    getOverlapMap() const;

    //! Get Jacobian graph
    virtual Teuchos::RCP<const Epetra_CrsGraph> 
    getJacobianGraph() const;

    //! Get overlap Jacobian graph
    virtual Teuchos::RCP<const Epetra_CrsGraph> 
    getOverlapJacobianGraph() const;

    //! Get number of nodes per element
    virtual int getNumNodesPerElement() const;


  private:

    //! Private to prohibit copying
    BlockDiscretization(const BlockDiscretization&);

    //! Private to prohibit copying
    BlockDiscretization& operator=(const BlockDiscretization&);

  protected:
    
    //! Underlying discretization object
    Teuchos::RCP<const FEApp::AbstractDiscretization> underlyingDisc;

    //! Epetra communicator
    Teuchos::RCP<const EpetraExt::MultiMpiComm> globalComm;

    //! Unknown Map
    Teuchos::RCP<Epetra_Map> map;

    //! Overlapped unknown map
    Teuchos::RCP<Epetra_Map> overlap_map;

    //! Jacobian matrix graph
    Teuchos::RCP<Epetra_CrsGraph> graph;

    //! Overlapped Jacobian matrix graph
    Teuchos::RCP<Epetra_CrsGraph> overlap_graph;

  };

}

#endif // FEAPP_CZERODISCRETIZATION_HPP
