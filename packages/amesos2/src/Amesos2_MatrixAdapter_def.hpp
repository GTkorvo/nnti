#ifndef AMESOS2_MATRIXADAPTER_DEF_HPP
#define AMESOS2_MATRIXADAPTER_DEF_HPP

namespace Amesos {

  
  template < class Matrix >
  MatrixAdapter<Matrix>::MatrixAdapter(RCP<Matrix> m)
    : mat_(m)
  {
    comm_ = static_cast<const adapter_t*>(this)->getComm_impl();
    col_map_ = static_cast<const adapter_t*>(this)->getColMap_impl();
    row_map_ = static_cast<const adapter_t*>(this)->getRowMap_impl();
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				const Teuchos::ArrayView<global_ordinal_t> colind,
				const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				typename MatrixAdapter<Matrix>::global_size_t& nnz,
				const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t, global_ordinal_t, node_t> > rowmap,
				Util::EStorage_Ordering ordering) const
  {
    help_getCrs(nzval, colind, rowptr,
		nnz, rowmap, ordering,
		typename adapter_t::get_crs_spec());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				const Teuchos::ArrayView<global_ordinal_t> colind,
				const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				typename MatrixAdapter<Matrix>::global_size_t& nnz,
				Util::EDistribution distribution,
				Util::EStorage_Ordering ordering) const
  {
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap
      = Util::getDistributionMap<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(distribution,
											this->getGlobalNumRows(),
											this->getComm());
    this->getCrs(nzval, colind, rowptr, nnz, Teuchos::ptrInArg(*rowmap), ordering);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				const Teuchos::ArrayView<global_ordinal_t> rowind,
				const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				typename MatrixAdapter<Matrix>::global_size_t& nnz,
				const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t, global_ordinal_t, node_t> > colmap,
				Util::EStorage_Ordering ordering) const
  {
    help_getCcs(nzval, rowind, colptr,
		nnz, colmap, ordering,
		typename adapter_t::get_ccs_spec());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				const Teuchos::ArrayView<global_ordinal_t> rowind,
				const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				typename MatrixAdapter<Matrix>::global_size_t& nnz,
				Util::EDistribution distribution,
				Util::EStorage_Ordering ordering) const
  {
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap
      = Util::getDistributionMap<local_ordinal_t,global_ordinal_t,global_size_t,node_t>(distribution,
											this->getGlobalNumCols(),
											this->getComm());
    this->getCcs(nzval, rowind, colptr, nnz, Teuchos::ptrInArg(*colmap), ordering);
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNumRows() const
  {
    return row_map_->getMaxAllGlobalIndex() + 1;
  }

  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNumCols() const
  {
    return col_map_->getMaxAllGlobalIndex() + 1;
  }
  
  template < class Matrix >
  typename MatrixAdapter<Matrix>::global_size_t
  MatrixAdapter<Matrix>::getGlobalNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getGlobalNNZ_impl();
  }
  
  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalNumRows() const
  {
    return row_map_->getNodeNumElements(); // TODO: verify
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalNumCols() const
  {
    return col_map_->getNodeNumElements();
  }
  
  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getLocalNNZ_impl();
  }

  template < class Matrix >
  std::string
  MatrixAdapter<Matrix>::description() const
  {
    std::ostringstream oss;
    oss << "Amesos2::MatrixAdapter wrapping: ";
    oss << mat_->description();
    return oss.str();
  }
  
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::describe(Teuchos::FancyOStream &out,
				  const Teuchos::EVerbosityLevel verbLevel) const
  {}
  


  /******************************
   * Private method definitions *
   ******************************/

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> colind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
				     Util::EStorage_Ordering ordering,
				     has_special_impl hsi) const
  {
    static_cast<const adapter_t*>(this)->getCrs_spec(nzval, colind, rowptr,
						     nnz, rowmap, ordering);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> colind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
				     Util::EStorage_Ordering ordering,
				     no_special_impl nsi) const
  {
    do_getCrs(nzval, colind, rowptr,
	      nnz, rowmap, ordering,
	      typename adapter_t::major_access());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> colind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
				   Util::EStorage_Ordering ordering,
				   row_access ra) const
  {
    using Teuchos::RCP;
    using Teuchos::ArrayView;
    using Teuchos::OrdinalTraits;
    
    RCP<const type> get_mat;
    if( *rowmap == *this->row_map_ ){
      // No need to redistribute
      get_mat = rcp(this,false); // non-owning
    } else {
      get_mat = get(rowmap);
    }
    // RCP<const type> get_mat = get(rowmap);

    // technically, if all is well and good, then rowmap == rmap
    RCP<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rmap = get_mat->getRowMap();
    TEUCHOS_ASSERT( *rowmap == *rmap );

    ArrayView<const global_ordinal_t> node_elements = rmap->getNodeElementList();
    typename ArrayView<const global_ordinal_t>::iterator row_it, row_end;
    row_end = node_elements.end();

    size_t rowptr_ind = OrdinalTraits<size_t>::zero();
    global_ordinal_t rowInd = OrdinalTraits<global_ordinal_t>::zero();
    for( row_it = node_elements.begin(); row_it != row_end; ++row_it ){
      rowptr[rowptr_ind++] = rowInd;
      size_t rowNNZ = get_mat->getGlobalRowNNZ(*row_it);
      size_t nnzRet = OrdinalTraits<size_t>::zero();
      ArrayView<global_ordinal_t> colind_view = colind.view(rowInd,rowNNZ);
      ArrayView<scalar_t> nzval_view = nzval.view(rowInd,rowNNZ);
      
      get_mat->getGlobalRowCopy(*row_it, colind_view, nzval_view, nnzRet);

      // It was suggested that instead of sorting each row's indices
      // individually, that we instead do a double-transpose at the
      // end, which would also lead to the indices being sorted.
      if( ordering == Util::Sorted_Indices ){
	Tpetra::sort2(colind_view.begin(), colind_view.end(), nzval_view.begin());
      }
      
      TEST_FOR_EXCEPTION( rowNNZ != nnzRet,
			  std::runtime_error,
			  "Number of values returned different from "
                          "number of values reported");
      rowInd += rowNNZ;
    }
    rowptr[rowptr_ind] = nnz = rowInd;
  }

  // TODO: This may not work with distributed matrices.
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCrs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> colind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> rowptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > rowmap,
				   Util::EStorage_Ordering ordering,
				   col_access ca) const
  {
    using Teuchos::Array;
    // get the ccs and transpose

    Array<scalar_t> nzval_tmp(nzval.size(), 0);
    Array<global_ordinal_t> rowind(colind.size(), 0);
    Array<global_size_t> colptr(this->getGlobalNumCols() + 1);
    this->getCcs(nzval_tmp(), rowind(), colptr(), nnz, rowmap, ordering);
    
    Util::transpose(nzval_tmp(), rowind(), colptr(), nzval, colind, rowptr);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> rowind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
				     Util::EStorage_Ordering ordering,
				     has_special_impl hsi) const
  {
    static_cast<const adapter_t*>(this)->getCcs_spec(nzval, rowind, colptr,
						     nnz, colmap, ordering);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::help_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				     const Teuchos::ArrayView<global_ordinal_t> rowind,
				     const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				     typename MatrixAdapter<Matrix>::global_size_t& nnz,
				     const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
				     Util::EStorage_Ordering ordering,
				     no_special_impl nsi) const
  {
    do_getCcs(nzval, rowind, colptr,
	      nnz, colmap, ordering,
	      typename adapter_t::major_access());
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> rowind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
				   Util::EStorage_Ordering ordering,
				   row_access ra) const
  {
    using Teuchos::Array;
    // get the crs and transpose

    Array<scalar_t> nzval_tmp(nzval.size(), 0);
    Array<global_ordinal_t> colind(rowind.size(), 0);
    Array<global_size_t> rowptr(this->getGlobalNumRows() + 1);
    this->getCrs(nzval_tmp(), colind(), rowptr(), nnz, colmap, ordering);
    
    Util::transpose(nzval_tmp(), colind(), rowptr(), nzval, rowind, colptr);
  }

  template < class Matrix >
  void
  MatrixAdapter<Matrix>::do_getCcs(const Teuchos::ArrayView<scalar_t> nzval,
				   const Teuchos::ArrayView<global_ordinal_t> rowind,
				   const Teuchos::ArrayView<typename MatrixAdapter<Matrix>::global_size_t> colptr,
				   typename MatrixAdapter<Matrix>::global_size_t& nnz,
				   const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > colmap,
				   Util::EStorage_Ordering ordering,
				   col_access ca) const
  {
    using Teuchos::RCP;
    using Teuchos::ArrayView;
    using Teuchos::OrdinalTraits;
    
    RCP<const type> get_mat;
    if( *colmap == *this->col_map_ ){
      // No need to redistribute
      get_mat = rcp(this,false); // non-owning
    } else {
      get_mat = get(colmap);
    }

    // If all is well and good, then colmap == cmap
    RCP<const Tpetra::Map<scalar_t,local_ordinal_t,global_ordinal_t> > cmap = get_mat->getColMap();
    TEUCHOS_ASSERT( *colmap == *cmap );

    ArrayView<global_ordinal_t> node_elements = cmap->getNodeElementList();
    typename ArrayView<global_ordinal_t>::iterator col_it, col_end;
    col_end = node_elements.end();

    size_t colptr_ind = OrdinalTraits<size_t>::zero();
    global_ordinal_t colInd = OrdinalTraits<global_ordinal_t>::zero();
    for( col_it = node_elements.begin(); col_it != col_end; ++col_it ){
      colptr[colptr_ind++] = colInd;
      size_t colNNZ = getGlobalColNNZ(*col_it);
      size_t nnzRet = 0;
      ArrayView<global_ordinal_t> rowind_view = rowind.view(colInd,colNNZ);
      ArrayView<scalar_t> nzval_view = nzval.view(colInd,colNNZ);
      getGlobalColCopy(*col_it, rowind_view, nzval_view, nnzRet);
      
      // It was suggested that instead of sorting each row's indices
      // individually, that we instead do a double-transpose at the
      // end, which would also lead to the indices being sorted.
      if( ordering == Util::Sorted_Indices ){
	Tpetra::sort2(rowind_view.begin(), rowind_view.end(), nzval_view.begin());
      }
      
      TEST_FOR_EXCEPTION( colNNZ != nnzRet,
			  std::runtime_error,
			  "Number of values returned different from "
                          "number of values reported");
      colInd += colNNZ;
    }
    colptr[colptr_ind] = nnz = colInd;
  }

  
  // These will link to concrete implementations
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getGlobalRowCopy(global_ordinal_t row,
					  const Teuchos::ArrayView<global_ordinal_t>& indices,
					  const Teuchos::ArrayView<scalar_t>& vals,
					  size_t& nnz) const
  {
    static_cast<const adapter_t*>(this)->getGlobalRowCopy_impl(row, indices, vals, nnz);
  }
  
  template < class Matrix >
  void
  MatrixAdapter<Matrix>::getGlobalColCopy(global_ordinal_t col,
					  const Teuchos::ArrayView<global_ordinal_t>& indices,
					  const Teuchos::ArrayView<scalar_t>& vals,
					  size_t& nnz) const
  {
    static_cast<const adapter_t*>(this)->getGlobalColCopy_impl(col, indices, vals, nnz);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getMaxRowNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getMaxRowNNZ_impl();
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getMaxColNNZ() const
  {
    return static_cast<const adapter_t*>(this)->getMaxColNNZ_impl();
  }
    
  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getGlobalRowNNZ(global_ordinal_t row) const
  {
    return static_cast<const adapter_t*>(this)->getGlobalRowNNZ_impl(row);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalRowNNZ(local_ordinal_t row) const
  {
    return static_cast<const adapter_t*>(this)->getLocalRowNNZ_impl(row);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getGlobalColNNZ(global_ordinal_t col) const
  {
    return static_cast<const adapter_t*>(this)->getGlobalColNNZ_impl(col);
  }

  template < class Matrix >
  size_t
  MatrixAdapter<Matrix>::getLocalColNNZ(local_ordinal_t col) const
  {
    return static_cast<const adapter_t*>(this)->getLocalColNNZ_impl(col);
  }

  template < class Matrix >
  bool
  MatrixAdapter<Matrix>::isLocallyIndexed() const
  {
    return static_cast<const adapter_t*>(this)->isLocallyIndexed_impl();
  }
  
  template < class Matrix >
  bool
  MatrixAdapter<Matrix>::isGloballyIndexed() const
  {
    return static_cast<const adapter_t*>(this)->isGloballyIndexed_impl();
  }
  
  template < class Matrix >
  RCP<const MatrixAdapter<Matrix> >
  MatrixAdapter<Matrix>::get(const Teuchos::Ptr<const Tpetra::Map<local_ordinal_t,global_ordinal_t,node_t> > map) const
  {
    return static_cast<const adapter_t*>(this)->get_impl(map);
  }
  

} // end namespace Amesos

#endif	// AMESOS2_MATRIXADAPTER_DEF_HPP
