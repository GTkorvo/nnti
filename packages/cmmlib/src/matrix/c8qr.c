#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

#define Q_BLOCK_DEFAULT 16

/*

   Build block QR factorization.

     A = QR
*/

void cmm_c8_qr(c8 *a, int as, int am, int an,
               c8 *q, int qs, int qm, int qn, int bq,
               c8 *r, int rs, int rm, int rn)
  {
   c8 *w;
   c8 *z;
   int i;

   if(am != qm)
     cmm_fatal("in cmm_c8_qr, am %d != qm %d\n",am,qm);

   if(an != rn)
     cmm_fatal("in cmm_c8_qr, an %d != rn %d\n",an,rn);

   if(rm != qn)
     cmm_fatal("in cmm_c8_qr, rm %d != qn %d\n",rm,qn);

   if(rm != rn)
     cmm_fatal("in cmm_c8_qr, rm %d != rn %d\n",rm,rn);


   if(bq == 0)
     bq = Q_BLOCK_DEFAULT;

// Allocate temp storage for W, Z
   w = cmm_alloc(sizeof(c8)*am*bq);
   z = cmm_alloc(sizeof(c8)*an*bq);

// Copy A to Q if not in-place
   if(q != a)
     cmm_c8_m_copy(a, as, am, an, q, qs);

   for(i=0;i<qn;i += bq)
     {
      c8 *qi;
      c8 *ar;
      c8 *ri;
      int bi, ni, mi, lr;

      mi = qm-i;
      ni = qn-i;
      bi = MIN(ni, bq);
 
      qi = q + i + qs*i;
      ri = r + i + rs*i;

// Build the block of Householder vectors
      c8_build_house(qi, qs, mi, bi, ri, rs); 

      lr = qn-(i+bi); // remaining number of columns to update
      if(lr <= 0)
        break;

// Build the associated WY transformation (Q = I+WY). Note
// that Y is just the Householder vectors.
      c8_build_wy(qi, qs, mi, bi, w, mi);

// Apply (I+WY)' to remainder of matrix.
      ar = q + i + (i+bi)*qs;  // pointer to remainder of matrix

// Compute Z = W'A
//      printf("i %d Z = W'*A\n",i);
      cmm_c8_mm_hn(w, mi, mi, bi,
                   ar, qs, mi, lr,
                   z,  bi, bi, lr);

// A = A + WZ
//      printf("i %d A = A + Y*Z\n",i);
      cmm_c8_mm_a_nn(qi, qs, mi, bi,
                     z,  bi, bi, lr,
                     ar, qs, mi, lr );

// Update R. Note that junk is left in upper part of Q.
// Should I clear it out?
      if(lr > 0)
        cmm_c8_m_copy(ar, qs, bi, lr, r + i + rs*(i+bi), rs);
      }

// Free temp storage
   cmm_free(w);
   cmm_free(z);
  }

/*

   Apply block Q, where Q is constructed by QR factorization.

   B = QA

*/

void cmm_c8_qapply(c8 *q, int qs, int qm, int qn, int bq,
                   c8 *a, int as, int am, int an,
                   c8 *b, int bs, int bm, int bn)
  {
   int i;
   c8 *w;
   c8 *z;

   if(an != bn)
     cmm_fatal("in cmm_c8_qr, an %d != bn %d\n",an,bn);

   if(am != qn)
     cmm_fatal("in cmm_c8_qr, am %d != qn %d\n",am,qn);

   if(bm != qm)
     cmm_fatal("in cmm_c8_qr, bm %d != qm %d\n",bm,qm);

   if(bq == 0)
     bq = Q_BLOCK_DEFAULT;

// Allocate temp storage for W, Z
   w = cmm_alloc(sizeof(c8)*am*bq);
   z = cmm_alloc(sizeof(c8)*qm*bq);

   for(i=0;i<qn;i+=bq)
     {
      c8 *qi;
      int bi, ni, mi, lr;

      mi = qm-i;
      ni = qn-i;
      bi = MIN(ni, bq);
 
// Get the block of Householder vectors
      qi = q + i + qs*i;      

// Build the WY representation
      c8_build_wy(qi, qs, mi, bi, w, mi);

// Apply (I+WY) to b

      lr = qn-(i+bi); // remaining number of columns to update
      if(lr <= 0)
        break;
     }


// Free temp storage
   cmm_free(w);
   cmm_free(z);

  }

/*

   Apply block Q', where Q is constructed by QR factorization.

   B = Q'A

*/

void cmm_c8_qtapply(c8 *q, int qs, int qm, int qn, int bq,
                   c8 *a, int as, int am, int an,
                   c8 *b, int bs, int bm, int bn)
  {

   if(an != bn)
     cmm_fatal("in cmm_c8_qr, an %d != bn %d\n",an,bn);

   if(am != qm)
     cmm_fatal("in cmm_c8_qr, am %d != qm %d\n",am,qm);

   if(bm != qn)
     cmm_fatal("in cmm_c8_qr, bm %d != qn %d\n",bm,qn);


   if(bq == 0)
     bq = Q_BLOCK_DEFAULT;


  }


/*

   Apply block Q, where Q is constructed by QR factorization.

   B = AQ

*/

void cmm_c8_applyq(c8 *q, int qs, int qm, int qn, int bq,
                  c8 *a, int as, int am, int an,
                  c8 *b, int bs, int bm, int bn)
  {

   if(am != bm)
     cmm_fatal("in cmm_c8_qr, am %d != bm %d\n",am,bm);

   if(an != qm)
     cmm_fatal("in cmm_c8_qr, an %d != qm %d\n",an,qm);

   if(bn != qn)
     cmm_fatal("in cmm_c8_qr, bn %d != qn %d\n",bn,qn);



   if(bq == 0)
    bq = Q_BLOCK_DEFAULT;


  }


/*

   Apply block Q', where Q is constructed by QR factorization.

   B = AQ'

*/

void cmm_c8_applyqt(c8 *q, int qs, int qm, int qn, int bq,
                   c8 *a, int as, int am, int an,
                   c8 *b, int bs, int bm, int bn)
  {

   if(am != bm)
     cmm_fatal("in cmm_c8_qr, am %d != bm %d\n",am,bm);

   if(an != qn)
     cmm_fatal("in cmm_c8_qr, an %d != qn %d\n",an,qn);

   if(bn != qm)
     cmm_fatal("in cmm_c8_qr, bn %d != qm %d\n",bn,qm);

   if(bq == 0)
    bq = Q_BLOCK_DEFAULT;

  }
