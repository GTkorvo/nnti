/*
 * MueLu_PgPFactory.hpp
 *
 *  Created on: 23.09.2011
 *      Author: tobias
 */

#ifndef MUELU_PGPFACTORY_HPP_
#define MUELU_PGPFACTORY_HPP_

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_TpetraVector.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
    @class PgPFactory class.
    @brief Factory for building Petrov-Galerkin Smoothed Aggregation prolongators.
 */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class PgPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, PgPFactory<AA,BB,CC,DD, EE> &factory);

public:

  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
   */
  PgPFactory(RCP<PFactory> InitialPFact = Teuchos::null, RCP<SingleLevelFactoryBase> AFact = Teuchos::null)
  : initialPFact_(InitialPFact), AFact_(AFact),
    diagonalView_("current") {
    if(initialPFact_ == Teuchos::null)
    {
      // use tentative P factory as default
      initialPFact_ = rcp(new TentativePFactory());
    }
  }

  //! Destructor.
  virtual ~PgPFactory() {}

  //@}

  //! @name Set methods.
  //@{

  //! Change view of diagonal.
  void SetDiagonalView(std::string const& diagView) {
    diagonalView_ = diagView;
  }
  //@}

  //! @name Get methods.
  //@{

  //! Returns current view of diagonal.
  std::string GetDiagonalView() {
    return diagonalView_;
  }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    //coarseLevel.Request("P",initialPFact_.get());
    //fineLevel.Request("A",AFact_.get());
    fineLevel.DeclareInput("A",AFact_.get());
    coarseLevel.DeclareInput("P",initialPFact_.get());
  };
  //@}

  //! @name Build methods.
  //@{

  /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
   */
  void Build(Level& fineLevel, Level &coarseLevel) const {
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
    std::ostringstream buf; buf << coarseLevel.GetLevelID();
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP_"+buf.str()));
    timer->start(true);

    // Level Get
    RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_.get());
    RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_.get());

    /////////////////// switch from A to A^T in restriction mode (necessary as long as implicit transpose not working for Epetra)
    if(restrictionMode_)
      A = Utils2<Scalar,LocalOrdinal,GlobalOrdinal>::Transpose(A,true); // build transpose of A explicitely

    /////////////////// calculate D^{-1} A Ptent (needed for smoothing)
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> DinvAP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

    doFillComplete=true;
    optimizeStorage=false;
    Teuchos::ArrayRCP<Scalar> diag = Utils::GetMatrixDiagonal(A);
    Utils::MyOldScaleMatrix(DinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    /////////////////// calculate local damping factors omega
    Teuchos::ArrayRCP<Scalar> RowBasedOmega = ComputeRowBasedOmegas(A,Ptent,DinvAP0,diag);


    /////////////////// prolongator smoothing using local damping parameters omega
    RCP<Operator> P_smoothed = Teuchos::null;
    Utils::MyOldScaleMatrix(DinvAP0,RowBasedOmega,false,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    Utils::TwoMatrixAdd(Ptent, false, Teuchos::ScalarTraits<Scalar>::one(),
                                         DinvAP0, false, -Teuchos::ScalarTraits<Scalar>::one(),
                                         P_smoothed);
    P_smoothed->fillComplete(Ptent->getDomainMap(), Ptent->getRangeMap());

    //////////////////// store results in Level

    // Level Set
    if(!restrictionMode_)
    {
        // prolongation factory is in prolongation mode
        coarseLevel.Set("P", P_smoothed, this);
    }
    else
    {
        // prolongation factory is in restriction mode
        RCP<Operator> R = Utils2<SC,LO,GO>::Transpose(P_smoothed,true); // use Utils2 -> specialization for double
        coarseLevel.Set("R", R, this);
    }

  }

  void BuildP(Level &fineLevel, Level &coarseLevel) const {
    std::cout << "TODO: remove me" << std::endl;
  }

  //@}

  Teuchos::ArrayRCP<Scalar> ComputeRowBasedOmegas(const RCP<Operator>& A, const RCP<Operator>& Ptent, const RCP<Operator>& DinvAPtent,const Teuchos::ArrayRCP<Scalar>& diagA) const
  {
    RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > colbasedomegamap = BuildLocalReplicatedColMap(DinvAPtent);

    if(colbasedomegamap->isDistributed() == true) std::cout << "colbasedomegamap is distributed" << std::endl; //throw("colbasedomegamap is distributed globally?");
    if(colbasedomegamap->getMinAllGlobalIndex() != DinvAPtent->getDomainMap()->getMinAllGlobalIndex()) std::cout << "MinAllGID does not match" << std::endl; //throw("MinAllGID does not match");
    if(colbasedomegamap->getMaxAllGlobalIndex() != DinvAPtent->getDomainMap()->getMaxAllGlobalIndex()) std::cout << "MaxAllGID does not match" << std::endl; //throw("MaxAllGID does not match");
    if(colbasedomegamap->getGlobalNumElements() != DinvAPtent->getDomainMap()->getGlobalNumElements()) std::cout << "NumGlobalElements do not match" << std::endl; //throw("NumGlobalElements do not match");


#if 0 // MUEMAT mode (=paper)
    // Minimize with respect to the (A)' A norm.
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' (A' A) D^{-1} A P0)
    //   omega =   ------------------------------------------
    //             diag( P0' A' D^{-1}' ( A'  A) D^{-1} A P0)
    //
    // expensive, since we have to recalculate AP0 due to the lack of an explicit scaling routine for DinvAP0

    // calculate A * Ptent
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> AP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

    // compute A * D^{-1} * A * P0
    RCP<Operator> ADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);

    RCP<Teuchos::Array<Scalar> > Numerator = MultiplyAll(AP0, ADinvAP0, colbasedomegamap);
    RCP<Teuchos::Array<Scalar> > Denominator = MultiplySelfAll(ADinvAP0, colbasedomegamap);
#elif 0
    // ML mode 1 (cheapest)
    // Minimize with respect to L2 norm
    //                  diag( P0' D^{-1} A P0)
    //   omega =   -----------------------------
    //             diag( P0' A' D^{-1}' D^{-1} A P0)
    //

    RCP<Teuchos::Array<Scalar> > Numerator = MultiplyAll(Ptent, DinvAPtent, colbasedomegamap);
    RCP<Teuchos::Array<Scalar> > Denominator = MultiplySelfAll(DinvAPtent, colbasedomegamap);
#elif 1
    // ML mode 2
    // Minimize with respect to the (D^{-1} A)' D^{-1} A norm.
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //   omega =   --------------------------------------------------------
    //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //

    // compute D^{-1} * A * D^{-1} * A * P0
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);
    Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    RCP<Teuchos::Array<Scalar> > Numerator = MultiplyAll(DinvAPtent, DinvADinvAP0, colbasedomegamap);
    RCP<Teuchos::Array<Scalar> > Denominator = MultiplySelfAll(DinvADinvAP0, colbasedomegamap);
#elif 0
    // ML mode 3 (most expensive)
    //             diag( P0' ( A'D' + DA) D A P0)
    //   omega =   -----------------------------
    //             diag( P0'A'D' ( A'D' + DA) D A P0)
    //
    //             diag( DinvAP0'DinvAP0 + P0'DinvADinvAP0)
    //         =   -----------------------------
    //                2*diag( DinvADinvAP0'DinvAP0)
    //
    //

    // compute D^{-1} * A * D^{-1} * A * P0
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);
    Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    RCP<Teuchos::Array<Scalar> > Numerator = MultiplyAll(Ptent, DinvADinvAP0, colbasedomegamap);
    RCP<Teuchos::Array<Scalar> > Numerator2= MultiplySelfAll(DinvAPtent, colbasedomegamap);
    if(Numerator->size() != Numerator2->size()) std::cout << "Error: size of Numerator and Numerator2 different. Error" << std::endl;
    for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
      (*Numerator)[i] += (*Numerator2)[i];
    RCP<Teuchos::Array<Scalar> > Denominator = MultiplyAll(DinvAPtent,DinvADinvAP0, colbasedomegamap);
    for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
      (*Denominator)[i] *= 2.;
#endif

    /////////////////// DEBUG: check for zeros in denominator
    size_t zeros_in_denominator = 0;
    for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
    {
      if((*Denominator)[i] == Teuchos::ScalarTraits<Scalar>::zero()) zeros_in_denominator ++;
    }
    if(zeros_in_denominator>Teuchos::ScalarTraits<Scalar>::zero()) std::cout << "There are " << zeros_in_denominator<< " zeros in Denominator. very suspicious!" << std::endl;

    /////////////////// build ColBasedOmegas
    RCP<Teuchos::ArrayRCP<Scalar> > ColBasedOmegas = Teuchos::rcp(new Teuchos::ArrayRCP<Scalar>(Numerator->size(),Teuchos::ScalarTraits<Scalar>::zero()));
    for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
    {
      (*ColBasedOmegas)[i] = (*Numerator)[i]/(*Denominator)[i];
      if((*ColBasedOmegas)[i] < Teuchos::ScalarTraits<Scalar>::zero())
        (*ColBasedOmegas)[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }

    /////////////////// transform ColBasedOmegas to row based omegas (local ids)
    Teuchos::ArrayRCP< Scalar > RowBasedOmegas =
        TransformCol2RowBasedOmegas(ColBasedOmegas, colbasedomegamap, DinvAPtent);

    return RowBasedOmegas;
  }

  // This routine still has Epetra/Tpetra specific code
  Teuchos::ArrayRCP< Scalar > TransformCol2RowBasedOmegas(const RCP<Teuchos::ArrayRCP<Scalar> >& ColBasedOmegas,const RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >& colbasedomegamap, const RCP<const Operator>& Op) const
  {
    RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RowBasedOmegas = Teuchos::null;

    // build RowBasedOmegas
    if (Op->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      RCP<Xpetra::EpetraVector> eRowBasedOmegas =
                rcp(new Xpetra::EpetraVector(Op->getRowMap(),true));
      eRowBasedOmegas->putScalar(-666.0);

      for(size_t row=0; row<Op->getNodeNumRows(); row++)
      {
        size_t nnz = Op->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const GlobalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Op->getLocalRowView(row, indices, vals);

        if(Teuchos::as<size_t>(indices.size()) != nnz) std::cout << "number of nonzeros not equal to number of indices?" << std::endl;

        if(nnz == 0)
        {
          eRowBasedOmegas->getEpetra_Vector()->ReplaceMyValue(row,0,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else if(nnz == 1)
        {
          // dirichlet dofs
          eRowBasedOmegas->getEpetra_Vector()->ReplaceMyValue(row,0,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else
        {
          for(size_t j=0; j<nnz; j++)
          {
            GlobalOrdinal col_gid = indices[j]; // local id = global id
            LocalOrdinal col_omega_lid = colbasedomegamap->getLocalElement(col_gid);
            Scalar omega = (*ColBasedOmegas)[col_omega_lid];
            Epetra_Vector* eeRowBasedOmegas = eRowBasedOmegas->getEpetra_Vector();
            if((*eeRowBasedOmegas)[row] == -666.0)
              eeRowBasedOmegas->ReplaceMyValue(row,0,omega);
            else if (omega < (*eeRowBasedOmegas)[row])
              eeRowBasedOmegas->ReplaceMyValue(row,0,omega);
          }
        }
      }
      RowBasedOmegas = eRowBasedOmegas;

#else
throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }

    // Tpetra-specific code
    else if(Op->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tRowBasedOmegas =
                rcp(new Xpetra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Op->getRowMap(),true));
      tRowBasedOmegas->putScalar(-666.0);
      for(size_t row=0; row<Op->getNodeNumRows(); row++)
      {
        size_t nnz = Op->getNumEntriesInLocalRow(row);

        Teuchos::ArrayView<const GlobalOrdinal> indices;
        Teuchos::ArrayView<const Scalar> vals;
        Op->getLocalRowView(row, indices, vals);

        if(Teuchos::as<size_t>(indices.size()) != nnz) std::cout << "number of nonzeros not equal to number of indices?" << std::endl;

        if(nnz == 0)
        {
          tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else if(nnz == 1)
        {
          // dirichlet dofs
          tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,Teuchos::ScalarTraits<Scalar>::zero());
        }
        else
        {
          for(size_t j=0; j<nnz; j++)
          {
            GlobalOrdinal col_gid = indices[j]; // local id = global id
            LocalOrdinal col_omega_lid = colbasedomegamap->getLocalElement(col_gid);
            Scalar omega = (*ColBasedOmegas)[col_omega_lid];
            RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ttRowBasedOmegas = tRowBasedOmegas->getTpetra_Vector();
            Teuchos::ArrayRCP<const Scalar> localRowData = ttRowBasedOmegas->getData(0);
            if(localRowData[row] == -666.0)
              tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,omega);
            else if (omega < localRowData[row])
              tRowBasedOmegas->getTpetra_Vector()->replaceLocalValue(row,omega);
          }
        }
      }

      RowBasedOmegas = tRowBasedOmegas;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }


    // check for negative entries in RowBasedOmegas
    Teuchos::ArrayRCP< Scalar > RowBasedOmega_data = RowBasedOmegas->getDataNonConst(Teuchos::as<size_t>(0));
    for(size_t i=0; i<Teuchos::as<size_t>(RowBasedOmega_data.size()); i++)
    {
      if(RowBasedOmega_data[i] < Teuchos::ScalarTraits<Scalar>::zero()) RowBasedOmega_data[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }

    return RowBasedOmegas->getDataNonConst(Teuchos::as<size_t>(0));
  }

  RCP<Teuchos::Array<Scalar> > MultiplySelfAll(const RCP<Operator>& Op, const RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >& InnerProdMap) const
  {
    Teuchos::Array<Scalar> InnerProd_local(InnerProdMap->getMaxAllGlobalIndex()+1,Teuchos::ScalarTraits<Scalar>::zero());

    for(size_t n=0; n<Op->getNodeNumRows(); n++)
    {
      Teuchos::ArrayView<const LocalOrdinal> lindices;
      Teuchos::ArrayView<const Scalar> lvals;
      Op->getLocalRowView(n, lindices, lvals);

      for(size_t i=0; i<Teuchos::as<size_t>(lindices.size()); i++)
      {
        GlobalOrdinal gid = Op->getColMap()->getGlobalElement(lindices[i]);
        InnerProd_local[gid] += lvals[i]*lvals[i];
      }
    }

    RCP<Teuchos::Array<Scalar> > InnerProd = rcp(new Teuchos::Array<Scalar>(InnerProdMap->getMaxAllGlobalIndex()+1,Teuchos::ScalarTraits<Scalar>::zero()));

    /////////////// sum up all values to global
    Teuchos::reduceAll(*(InnerProdMap->getComm()),Teuchos::REDUCE_SUM, Teuchos::as<const GlobalOrdinal>(InnerProd->size()) ,&InnerProd_local[0], &(*InnerProd)[0]);

    return InnerProd;
  }



  RCP<Teuchos::Array<Scalar> > MultiplyAll(const RCP<Operator>& left, const RCP<Operator>& right, const RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >& InnerProdMap) const
  {

    if(!left->getDomainMap()->isSameAs(*right->getDomainMap()))
      std::cout << "domain maps of left and right do not match" << std::endl;
    if(!left->getRowMap()->isSameAs(*right->getRowMap()))
      std::cout << "row maps of left and right do not match" << std::endl;

    // test
    /*Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
    Teuchos::RCP<Xpetra::Operator<double> > leftT = MueLu::Utils2<double>::Transpose(left,false);
    Teuchos::RCP<Xpetra::Operator<double> > leftTright = MueLu::Utils<double>::TwoMatrixMultiply(leftT,false,right,false);
    leftTright->describe(*fos,Teuchos::VERB_EXTREME);*/

    Teuchos::Array<Scalar> InnerProd_local(InnerProdMap->getMaxAllGlobalIndex()+1,Teuchos::ScalarTraits<Scalar>::zero());

    for(size_t n=0; n<left->getNodeNumRows(); n++)
    {
      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      left->getLocalRowView(n, lindices_left, lvals_left);

      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;
      right->getLocalRowView(n, lindices_right, lvals_right);

      for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
      {
        for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
        {
          GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
          GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
          if(left_gid == right_gid)
          {
            InnerProd_local[left_gid] += lvals_left[i]*lvals_right[j];
            break; // skip remaining gids of right operator
          }
        }
      }
    }

    RCP<Teuchos::Array<Scalar> > InnerProd = rcp(new Teuchos::Array<Scalar>(InnerProdMap->getMaxAllGlobalIndex()+1,Teuchos::ScalarTraits<Scalar>::zero()));

    /////////////// sum up all values to global
    Teuchos::reduceAll(*(InnerProdMap->getComm()),Teuchos::REDUCE_SUM, Teuchos::as<const GlobalOrdinal>(InnerProd->size()) ,&InnerProd_local[0], &(*InnerProd)[0]);

    return InnerProd;
  }



  Teuchos::RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > BuildLocalReplicatedColMap(const RCP<Operator>& DinvAP0) const
  {
    //Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(cout));
    Teuchos::RCP< const Teuchos::Comm< int > > comm = DinvAP0->getRangeMap()->getComm();

    // generate allreduced column-based array vector for local damping factors
    Teuchos::Array<GlobalOrdinal> localreplicatedcolgids;
    reduceAllXpetraMap(localreplicatedcolgids,*(DinvAP0->getDomainMap()));

    // Epetra-specific code
    if (DinvAP0->getRowMap()->lib() == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
      Teuchos::RCP<Epetra_Map> epetra_locreplicatedomegagids = Teuchos::rcp(new Epetra_Map(localreplicatedcolgids.size(),localreplicatedcolgids.size(),&localreplicatedcolgids[0],0,*Xpetra::toEpetra(comm)));
      const RCP< const Xpetra::EpetraMap> xpetra_locreplicatedomegagids = Teuchos::rcp(new Xpetra::EpetraMap(epetra_locreplicatedomegagids));
      return xpetra_locreplicatedomegagids;
#else
      throw("HAVE_MUELU_EPETRA_AND_EPETRAEXT not set. Compile MueLu with Epetra and EpetraExt enabled.")
#endif
    }
    // Tpetra-specific code
    else if(DinvAP0->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      Teuchos::RCP<Teuchos::ArrayView<const GlobalOrdinal> > arView = Teuchos::rcp(new Teuchos::ArrayView<const GlobalOrdinal>(&localreplicatedcolgids[0],localreplicatedcolgids.size()));
      //Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tpetra_locreplicatedomegagids = Tpetra::createNonContigMap<GlobalOrdinal,LocalOrdinal>(*arView, comm);
      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tpetra_locreplicatedomegagids = Tpetra::createLocalMap<GlobalOrdinal,LocalOrdinal>(Teuchos::as<size_t>(DinvAP0->getColMap()->getMaxAllGlobalIndex()+1), comm); // create local map
      const RCP< const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > xpetra_locreplicatedomegagids = Xpetra::toXpetra(tpetra_locreplicatedomegagids);
      return xpetra_locreplicatedomegagids;
#else
      throw("HAVE_MUELU_TPETRA not set. Compile MueLu with Tpetra enabled.")
#endif
    }

    return Teuchos::null;
  }


private:

  //! @name helper function for allreducing a Xpetra::Map
  //@{

  int FindMyPos(size_t nummyelements, const Teuchos::Comm<int>& comm) const
  {
    std::cout << "begin FindMyPos" << std::endl;
    const int myrank = comm.getRank();
    const int numproc = comm.getSize();

    std::vector<size_t> snum (numproc,0); // local vector of no of elements
    std::vector<size_t> rnum (numproc);   // gobal vector of no of elements on proc
    snum[myrank] = nummyelements;

    std::cout << "FindMyPose::reduce all" << std::endl;
    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, numproc, &snum[0], &rnum[0]);
    std::cout << "FindMyPos:: after reduce all" << std::endl;

    return std::accumulate(&rnum[0], &rnum[myrank], 0);
  }

  void reduceAllXpetraMap(Teuchos::Array<GlobalOrdinal>& rredundant, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const
  {
    const int mynodepos = FindMyPos(map.getNodeNumElements(), *(map.getComm()));

    std::vector<GlobalOrdinal> sredundant(map.getGlobalNumElements(),0);

    Teuchos::ArrayView< const GlobalOrdinal > gids = map.getNodeElementList();
    std::copy(gids.getRawPtr(), gids.getRawPtr() + map.getNodeNumElements(), &sredundant[mynodepos]); // check me

    rredundant.resize(map.getGlobalNumElements());
    GlobalOrdinal numEle = map.getGlobalNumElements();
    std::cout << "reduce All in XpetraMap" << std::endl;
    Teuchos::reduceAll(*(map.getComm()),Teuchos::REDUCE_SUM, numEle ,&sredundant[0], &rredundant[0]);
    std::cout << "after reduce All in XpetraMap" << std::endl;
  }

  //@}


private:

  //! Input factories
  RCP<PFactory> initialPFact_;        //! Ptentative Factory
  RCP<SingleLevelFactoryBase> AFact_; //! A Factory

  //! Factory parameters
  std::string diagonalView_;
};

//! Friend print function.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, PgPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
  os << "Printing an PgPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_PGPFACTORY_SHORT

#endif /* MUELU_PGPFACTORY_HPP_ */
