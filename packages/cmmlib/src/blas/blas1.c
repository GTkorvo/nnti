/*! \file blas1.c
   \brief Level 1 BLAS and miscellaneous routines
   \ingroup blas
*/

#include "cmm.h"
#include "cmmblas.h"

/*! \brief Clear integer vector.

    \param n length of vector.
    \param iv integer vector.

*/
void cmm_i_zero(int n, int *iv)
    {
     int i;

     for(i=0;i<n;i++)
        *iv++ = 0;
    }

/*! \brief Set an integer vector to a constant.

    \param n length of vector.
    \param c constant value.
    \param iv integer vector.
*/
void cmm_i_set(int n, int c, int *iv)
    {
     int i;

     for(i=0;i<n;i++)
        *iv++ = c;
    }

/*! \brief Copy vector to another vector with strides.

    \param n Number of elements to copy.
    \param src Source vector.
    \param src_stride Source stride.
    \param dst Destination vector.
    \param dst_stride Destination stride.

*/
void cmm_i_copy(int n, int *src, int src_stride, int *dst, int dst_stride)
    {
     int i;

     for(i=0;i<n;i++)
       {
        *dst = *src;
         src += src_stride;
         dst += dst_stride;
       }
    }

/*! \brief Gather integer vector elements to a packed vector.

    \param n Number of elements in list.
    \param list List of elements to gather.
    \param src  Source vector.
    \param dst Destination (packed) vector.

*/
void cmm_i_gather(int n, int *list, int *src, r8 *dst)
    {
     int i;

     for(i=0;i<n;i++)
       {
        int j = *list++;
        *dst++ = src[j];
       }
    }

/*! \brief Gather integer vector blocks to packed vector.

    \param n    Number of blocks.
    \param block Number of elements per block.
    \param list  List of blocks to gather.
    \param src  Source vector.
    \param dst Destination vector.
*/
void cmm_i_gather_block(int n, int block, int *list, int *src, r8 *dst)
    {
     int i;

     for(i=0;i<n;i++)
       {
        int j = *list++;
        int k;
        int *ju = src+block*j;
        for(k=0;k<block;k++)
          *dst++ = *ju++;
       }
    }

/*! \brief Scallter elements of integer vector.

    \param n Number of elements to scatter.
    \param list List of elements to scatter.
    \param src Source (packed) vector.
    \param dst Destination vector.

*/
void cmm_i_scatter(int n, int *list, int *src, r8 *dst)
    {
     int i;

     for(i=0;i<n;i++)
       {
        int j = *list++;
        dst[j] = *src++;
       }
    }

/*! \brief Find minimum and maximum elements of an integer vector.

    \param n Number of elements.
    \param iv Vector.
    \param ivmin  Minimum value.
    \param imin   Location of minimum value.
    \param ivmax  Maximum value.
    \param imax   Location of maximum value.

*/
void cmm_i_range(int n, int *iv, int *ivmin, int *imin, int *ivmax, int *imax)
    {
     int i;
     int v0, v1, i0, i1;

     i0 = 0;
     i1 = 0;

     v0 = iv[0];
     v1 = iv[0];
     iv++;

     for(i=1;i<n;i++)
       {
        int v = *iv++;

        if(v < v0)
          {
           v0 = v;
           i0 = i;
          }

        if(v > v1)
          {
           v1 = v;
           i1 = i;
          }

       }

     if(ivmin)
       *ivmin = v0;
     if(imin)
       *imin  = i0;

     if(ivmax)
       *ivmax = v1;
     if(imax)
       *imax  = i1;
    }

/*! \brief Sort integer vector.

    \param n Number of elements.
    \param iv Vector.
    \param list Sort order.

*/
void cmm_i_sort(int n, int *iv, int *list)
    {
     cmm_fatal("cmm_i_sort not implemented\n");
    }

// REAL ROUTINES

/*! \brief Clear real vector.

    \param n Number of elements.
    \param v Vector to clear.

*/
void cmm_r8_zero(int n, r8 *v)
    {
     int i;

     for(i=0;i<n;i++)
        *v++ = 0.;
    }

/*! \brief Set real vector to constant.

    \param n Number of elements.
    \param c Value to set.
    \param v Vector.

*/
void cmm_r8_set(int n, r8 c, r8 *v)
    {
     int i;

     for(i=0;i<n;i++)
        *v++ = c;
    }

/*! \brief Copy a real vector to another vector, with strides.

    \param n Number of elements to copy.
    \param src Source vector.
    \param src_stride Source stride.
    \param dst Destination vector.
    \param dst_stride Destination stride.

*/
void cmm_r8_copy(int n, r8 *src, int src_stride, r8 *dst, int dst_stride)
    {
     int i;

     for(i=0;i<n;i++)
       {
        *dst = *src;
         src += src_stride;
         dst += dst_stride;
       }
    }

/*! \brief Gather elements  of a vector.

    \param n         number of elements.
    \param list      list of elements to gather
    \param src       source (unpacked) vector
    \param dst       destination (packed) vector

*/
void cmm_r8_gather(int n, int *list, r8 *src, r8 *dst)
    {
     int i;

     for(i=0;i<n;i++)
       {
        int j = *list++;
        *dst++ = src[j];
       }
    }

/*! \brief Scatter elements of a vector.

    \param n          number of elements.
    \param list      list of elements to scatter
    \param src       source (packed) vector
    \param dst       destination (unpacked) vector

*/
void cmm_r8_scatter(int n, int *list, r8 *src, r8 *dst)
    {
     int i;

     for(i=0;i<n;i++)
       {
        int j = *list++;
        dst[j] = *src++;
       }
    }

/*! \brief sort a real vector.

    \param n        number of elements.
    \param v        vector
    \param list     sort order

*/
void cmm_r8_sort(int n, r8 *v, int *list)
    {
     cmm_fatal("cmm_r8_sort not implemented\n");
    }

/*!
    \brief scale a vector

    \param n          number of elements
    \param   c          scale factor
    \param v          vector to scale

*/
void cmm_r8_scale(int n, r8 c, r8 *v)
    {
     int i;

     for(i=0;i<n;i++)
       *v++ *= c;

    }

/*!    \brief compute vector sum

    \param n          number of elements
    \param v          vector to sum

    \return Sum of vector elements.

*/
r8   cmm_r8_sum(int n, r8 *v)
    {
     int i;
     r8 t;

     t = 0.;
     for(i=0;i<n;i++)
       t += *v++;

     return t;
    }

/*!    \brief compute vector dot product

    \param n          number of elements
    \param u          vector
    \param v          vector

    \return dot product of vectors.

*/
r8   cmm_r8_dot(int n, r8 *u, r8 *v)
    {
     int i;
     r8 t;

     t = 0.;
     for(i=0;i<n;i++)
       t += *u++ * *v++;

     return t;
    }

/*!    \brief compute vector norm

    \param n          number of elements
    \param u          vector

    \return norm of vector

*/
r8   cmm_r8_norm(int n, r8 *u)
    {
     return sqrt(cmm_r8_normsq(n,u));
    }

/*!    \brief compute vector norm squared

    \param n          number of elements
    \param u          vector

    \return norm squared of vector

*/
r8   cmm_r8_normsq(int n, r8 *u)
    {
     int i;
     r8 t;

     t = 0.;
     for(i=0;i<n;i++)
       {
        r8 x;
        x = *u++;
        t += x*x;
       }

     return t;
    }

/*!
    \brief compute vector sum of absolute values (1-norm)

    \param n          number of elements
    \param u          vector

    \return 1-norm of vector

*/
r8   cmm_r8_norm1(int n, r8 *u)
    {
     int i;
     r8 t;

     t = 0.;
     for(i=0;i<n;i++)
       {
        r8 x;
        x = *u++;
        t += ABS(x);
       }

     return t;
    }

/*!
    \brief compute infinity norm of vector

    \param n          number of elements
    \param u          vector

    \return infinity norm of vector (max absolute value)

*/
r8   cmm_r8_normi(int n, r8 *u)
    {
     int i;
     r8 t;

     t = ABS(u[0]);
     for(i=0;i<n;i++)
       {
        r8 x;
        x = *u++; x = ABS(x);  // Danger, will robinson! be careful with ABS.
        t = MAX(t, x);
       }

     return t;
    }

/*!    \brief find min, max elements of a vector

    \param n          number of elements.
    \param v         vector
    \param vmin      minimum value returned
    \param imin      minimum value location
    \param vmax      maximum value returned
    \param imax      maximum value location

*/
void cmm_r8_range(int n, r8 *v, r8 *vmin, int *imin, r8 *vmax, int *imax)
    {
     int i;
     int i0, i1;
     r8 v0, v1;

     i0 = 0;
     i1 = 0;

     v0 = v[0];
     v1 = v[0];
     v++;

     for(i=1;i<n;i++)
       {
        r8 t = *v++;

        if(t < v0)
          {
           v0 = t;
           i0 = i;
          }

        if(t > v1)
          {
           v1 = t;
           i1 = i;
          }

       }

     if(vmin)
       *vmin = v0;
     if(imin)
       *imin  = i0;

     if(vmax)
       *vmax = v1;
     if(imax)
       *imax  = i1;
    }


/*!    \brief y = y + a*x (AXPY operation).

    \param n     number of elements.
    \param a      scalar factor
    \param x     vector to add
    \param y     vector result

*/
void cmm_r8_axpy(int n, r8 a, r8 *x, r8 *y)
    {
     int i;

     for(i=0;i<n;i++)
       *y++ += a * *x++;
    }

/*!
    \brief copy vector elements.

    \param n
    \param src
    \param src_stride
    \param dst
    \param dst_stride

*/

void cmm_c8_copy(int n, c8 *src, int src_stride, c8 *dst, int dst_stride)
    {
     int i;

     for(i=0;i<n;i++)
       {
        *dst = *src;
         src += src_stride;
         dst += dst_stride;
       }
    }

/*!
    \brief gather vector elements

    \param n        number of elements to gather
    \param list     list of elements
    \param u        source vector     
    \param v        destination (packed) vector

    
*/

void cmm_c8_gather(int n, int *list, c8 *u, c8 *v)
    {
     int i;

     for(i=0;i<n;i++)
       {
         int j = *list++;
        *v++ = u[j];
       }
    }

/*!    \brief scatter vector elements

    \param n        Number of elements to scatter.
    \param list     List of elements.
    \param u        Source (packed) vector.
    \param v        Destination vector.

*/

void cmm_c8_scatter(int n, int *list, c8 *u, c8 *v)
    {
     int i;

     for(i=0;i<n;i++)
       {
         int j = *list++;
         v[j] = *u++;
       }
    }

/*!    \brief clear complex vector.
       \param n Number of elements.
       \param v Vector.
    
*/

void cmm_c8_zero(int n, c8 *v)
    {
     int i;

     for(i=0;i<n;i++)
       {
        v->r = 0.;
        v->i = 0.;
        v++;
       }

    }


/*!
    \brief set complex vector to constant.
    \param n Number of elements.
    \param c Constant value.
    \param v Vector.
    
*/
void cmm_c8_set(int n, c8 c, c8 *v)
    {
     int i;
     for(i=0;i<n;i++)
       v[i] = c;
    }

/*!    \brief scale complex vector.
      \param n Number of elements.
    \param c Scale factor.
    \param v Vector.
    
*/
void cmm_c8_scale(int n, c8 c, c8 *v)
    {
     int i;
     r8 ar;
     r8 ai;

     for(i=0;i<n;i++)
       {
        ar = c.r*v->r - c.i*v->i;
        ai = c.r*v->i + c.i*v->r;
        v->r = ar;
        v->i = ai;
        v++;
       }
    }

/*!    \brief sum complex vector.

     \param n Number of elements.
     \param v Vector.
    \return sum of elements.

*/
c8   cmm_c8_sum(int n, c8 *v)
    {
     int i;
     c8 z;
     r8 ar = 0.;
     r8 ai = 0.;

     for(i=0;i<n;i++)
       {
        ar += v->r;
        ai += v->i;
        v++;
       }

     z.r = ar;
     z.i = ai;
     return z;
    }

/*!    \brief dot product of two complex vectors.

   \param n Number of elements.
   \param x First vector.
   \param y Second vector.
    \return complex dot product.

*/

c8   cmm_c8_dot(int n, c8 *x, c8 *y)
    {
     int i;
     c8 z;
     r8 ar = 0.;
     r8 ai = 0.;

     for(i=0;i<n;i++)
       {
        ar += x->r*y->r - x->i*y->i;
        ai += x->r*y->i + x->i*y->r;
        x++; y++;
       }

     z.r = ar;
     z.i = ai;
     return z;
    }

/*!
    \brief hermitian (conjugate) dot product

   \param n Number of elements.
   \param x First vector.
   \param y Second vector.
    \return hermitian dot product. Note that this is x dot conj(y)

*/

c8   cmm_c8_hdot(int n, c8 *x, c8 *y)
    {
     int i;
     c8 z;
     r8 ar = 0.;
     r8 ai = 0.;

     for(i=0;i<n;i++)
       {
        ar += x->r*y->r + x->i*y->i;
        ai += x->i*y->r - x->r*y->i;
        x++; y++;
       }

     z.r = ar;
     z.i = ai;
     return z;
    }


/*!
    \brief compute norm of complex vector.

   \param n Number of elements.
   \param v Vector.
    \return  value of norm.

*/

r8   cmm_c8_norm(int n, c8 *v)
    {
     return sqrt(cmm_c8_normsq(n,v));
    }

/*!
    \brief compute norm squared of complex vector.

   \param n Number of elements.
   \param x Vector.

    \return value of norm squared.

*/

r8   cmm_c8_normsq(int n, c8 *x)
    {
     int i;
     r8 t;

     t = 0.;

     for(i=0;i<n;i++)
       {
        t += x->r*x->r + x->i*x->i;
        x++;
       }

     return t;
    }


/*!
    \brief a*x+y

    \param n     number of elements.
    \param  a      scalar factor
    \param x     vector to add
    \param y     vector result

*/

void cmm_c8_axpy(int n, c8 a, c8 *x, c8 *y)
    {
     int i;
     r8 ar = a.r;
     r8 ai = a.i;

     for(i=0;i<n;i++)
       {
        y->r += ar*x->r - ai*x->i;
        y->i += ai*x->r + ar*x->i;
        x++; y++;
       }
    }

