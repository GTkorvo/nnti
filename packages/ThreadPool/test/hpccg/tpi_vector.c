
#include <stddef.h>

#include <TPI.h>
#include <tpi_vector.h>

/*--------------------------------------------------------------------*/

struct tpi_work_vector {
        double alpha ;
        double beta ;
  const double * x ;
  const double * y ;
        double * w ; 
        int  n ;
};

static void tpi_work_span( TPI_Work * const work , const int n ,
                           int * const iBeg , int * const iEnd )
{
  const int chunk = ( n + work->count - 1 ) / work->count ;

  *iEnd = chunk * ( work->rank + 1 );
  *iBeg = chunk * ( work->rank );

  if ( n < *iEnd ) { *iEnd = n ; }
}

/*--------------------------------------------------------------------*/

static void tpi_work_fill( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const double alpha = h->alpha ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = alpha ; }
}

void tpi_fill( int n , double alpha , double * x )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.w = x ;
  tmp.n = n ;
  TPI_Run_threads( tpi_work_fill , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_scale( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const double beta = h->beta ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] *= beta ; }
}

void tpi_scale( int n , const double alpha , double * x )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.w = x ;
  tmp.n = n ;
  TPI_Run_threads( tpi_work_scale , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_copy( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const double * const x = h->x ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = x[i] ; }
}

void tpi_copy( int n , const double * x , double * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.x = x ;
  tmp.w = y ;
  tmp.n = n ;
  TPI_Run_threads( tpi_work_copy , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_axpby( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const double alpha = h->alpha ;
  const double beta  = h->beta ;
  const double * const x = h->x ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = alpha * x[i] + beta * w[i] ; }
}

void tpi_axpby( int n , double alpha , const double * x ,
                        double beta  ,       double * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.beta  = beta ;
  tmp.x = x ;
  tmp.w = y ;

  if ( 0.0 == alpha ) {
    TPI_Run_threads( tpi_work_scale , & tmp , 0 );
  }
  else if ( 0.0 == beta && 1.0 == alpha ) {
    TPI_Run_threads( tpi_work_copy , & tmp , 0 );
  }
  else {
    TPI_Run_threads( tpi_work_axpby , & tmp , 0 );
  }
}

/*--------------------------------------------------------------------*/

static void tpi_work_dot_partial( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  double * const s = (double *) work->reduce ;
  const double * const x = h->x ;
  const double * const y = h->y ;
  double tmp = *s ;
  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { tmp += x[i] * y[i] ; }

  *s = tmp ;
}

static void tpi_work_dot_partial_self( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  double * const s = (double *) work->reduce ;
  const double * const x = h->x ;
  double tmp = *s ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { const double d = x[i] ; tmp += d * d ; }

  *s = tmp ;
}

static void tpi_work_dot_reduce( TPI_Work * work , const void * src  )
{
  *((double *) ( work->reduce) ) += *((const double *) src);
}

double tpi_dot( int n , const double * x , const double * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  double result = 0.0 ;
  tmp.x = x ;
  tmp.y = y ;
  tmp.n = n ;
  if ( x != y ) {
    TPI_Run_threads_reduce( tpi_work_dot_partial , & tmp ,
                            tpi_work_dot_reduce , & result , sizeof(result) );
  }
  else {
    TPI_Run_threads_reduce( tpi_work_dot_partial_self , & tmp ,
                            tpi_work_dot_reduce , & result , sizeof(result) );
  }
  return result ;
}

/*--------------------------------------------------------------------*/

struct tpi_crs_matrix {
        int      nRow ;
  const int    * A_pc ;
  const int    * A_ia ;
  const float  * A_a ;
  const double * x ;
        double * y ;
};

static void tpi_work_crs_matrix_apply( TPI_Work * work )
{
  const struct tpi_crs_matrix * const h =
    (struct tpi_crs_matrix *) work->info ;

  const int   * const A_pc = h->A_pc ;
  const int   * const A_ia = h->A_ia ;
  const float * const A_a  = h->A_a ;
  const double * const x = h->x ;
        double * const y = h->y ;

  const int nRow  = h->nRow ;
  const int chunk = ( nRow + work->count - 1 ) / work->count ;

  int rowEnd = chunk * ( work->rank + 1 );
  int row    = chunk * work->rank ;

  if ( nRow < rowEnd ) { rowEnd = nRow ; }

  for ( ; row < rowEnd ; ++row ) {
    const int jEnd = A_pc[ row + 1 ];
    int j = A_pc[ row ];
    double tmp = 0 ;
    for ( ; j < jEnd ; ++j ) { tmp += A_a[j] * x[ A_ia[j] ]; }
    y[ row ] = tmp ;
  }
}

/*--------------------------------------------------------------------*/

void tpi_crs_matrix_apply(
  const int      nRow ,
  const int    * A_pc ,
  const int    * A_ia ,
  const float  * A_a ,
  const double * x ,
        double * y )
{
  struct tpi_crs_matrix h = { 0 , NULL , NULL , NULL , NULL , NULL };
  h.nRow = nRow ;
  h.A_pc = A_pc ;
  h.A_ia = A_ia ;
  h.A_a  = A_a ;
  h.x    = x ;
  h.y    = y ;
  TPI_Run_threads( tpi_work_crs_matrix_apply , & h , 0 );
}


