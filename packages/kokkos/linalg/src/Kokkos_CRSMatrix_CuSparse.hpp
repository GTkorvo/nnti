#include <impl/Kokkos_PhysicalLayout.hpp>
namespace Kokkos {


/// CuSparse only works for double and float as scalar_type: return for all other scalar_types
template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_DoCuSparse(typename Kokkos::Impl::enable_if<
		!Kokkos::Impl::is_same<T,double>::value && !Kokkos::Impl::is_same<T,float>::value, typename RangeVector::scalar_type  >::type s_b
		,const RangeVector & y, typename DomainVector::scalar_type s_a,
		const CrsMatrix & A , const DomainVector & x) {
	return false;
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_DoCuSparse(typename Kokkos::Impl::enable_if<Kokkos::Impl::is_same<T,double>::value, double  >::type s_b
		,const RangeVector & y, double s_a,
		const CrsMatrix & A , const DomainVector & x) {

  /// If input is simple vectors, don't check for Layout Compatibility
  if(x.dimension_1()==1)
  cusparseDcsrmv(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                 A.numRows(), A.numCols(),  A.nnz(),
                 &s_a,
                 A.cusparse_descr,
                 A.values.ptr_on_device(),
                 (const int*) A.graph.row_map.ptr_on_device(),
                 A.graph.entries.ptr_on_device(),
                 x.ptr_on_device(),
                 &s_b,
                 y.ptr_on_device());
  else {
 
    /// Check that layouts of multivectors are compatible to what CuSparse expects
    Impl::PhysicalLayout layout_x(x);
    Impl::PhysicalLayout layout_y(y);
    if((layout_x.layout_type!=layout_x.Left) || layout_y.layout_type!=layout_y.Left) return false;
	
    cusparseDcsrmm(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                   A.numRows(), x.dimension_1(), A.numCols(),  A.nnz(),
                   &s_a,
                   A.cusparse_descr,
                   A.values.ptr_on_device(),
                   (const int*) A.graph.row_map.ptr_on_device(),
                   A.graph.entries.ptr_on_device(),
                   x.ptr_on_device(),
                   layout_x.stride[1],
                   &s_b,
                   y.ptr_on_device(),
                   layout_y.stride[1]);
  }
  return true;
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_DoCuSparse(typename Kokkos::Impl::enable_if<Kokkos::Impl::is_same<T,float>::value, float  >::type s_b
		,const RangeVector & y, float s_a,
		const CrsMatrix & A , const DomainVector & x) {
	
  if(x.dimension_1()==1)
    cusparseScsrmv(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                   A.numRows(), A.numCols(),  A.nnz(),
                   &s_a,
                   A.cusparse_descr,
                   A.values.ptr_on_device(),
                   (const int*) A.graph.row_map.ptr_on_device(),
                   A.graph.entries.ptr_on_device(),
                   x.ptr_on_device(),
                   &s_b,
                   y.ptr_on_device());
  else {
    
    /// Check that layouts of multivectors are compatible to what CuSparse expects
    Impl::PhysicalLayout layout_x(x);
    Impl::PhysicalLayout layout_y(y);
    if((layout_x.layout_type!=layout_x.Left) || layout_y.layout_type!=layout_y.Left) return false;

    cusparseScsrmm(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                   A.numRows(), x.dimension_1(), A.numCols(),  A.nnz(),
                   &s_a,
                   A.cusparse_descr,
                   A.values.ptr_on_device(),
                   (const int*) A.graph.row_map.ptr_on_device(),
                   A.graph.entries.ptr_on_device(),
                   x.ptr_on_device(),
                   layout_x.stride[1],
                   &s_b,
                   y.ptr_on_device(),
                   layout_y.stride[1]);	
  }
  return true;
}

//ToDo: strip compatible type attributes (const, volatile); make type of s_b and s_a independent
template<class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_Try_CuSparse( typename RangeVector::scalar_type s_b,const RangeVector & y, typename DomainVector::scalar_type s_a,
		const CrsMatrix & A , const DomainVector & x)
{
  if(!Kokkos::Impl::is_same<typename RangeVector::device_type,typename Kokkos::Cuda>::value) return false;
  if(Kokkos::Impl::is_same<typename RangeVector::non_const_scalar_type,float>::value&&
	 Kokkos::Impl::is_same<typename DomainVector::non_const_scalar_type,float>::value&&
	 Kokkos::Impl::is_same<typename CrsMatrix::values_type::non_const_scalar_type,float>::value) {
	   return MV_Multiply_DoCuSparse<typename RangeVector::scalar_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
  } else
  if(Kokkos::Impl::is_same<typename RangeVector::non_const_scalar_type,double>::value&&
	 Kokkos::Impl::is_same<typename DomainVector::non_const_scalar_type,double>::value&&
	 Kokkos::Impl::is_same<typename CrsMatrix::values_type::non_const_scalar_type,double>::value) {
	   return MV_Multiply_DoCuSparse<typename RangeVector::scalar_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
  } else
  return false;
}

}
