      subroutine scsrmm2( m, n, val, indx, pntr, x, y, nrhs)
*
*     Performs the matrix-vector operation
*
*                               y = A*x
*
*     where x and y are vectors and A is a sparse matrix stored
*     in compress row format.
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      implicit none
      integer m, n, nrhs, indx(*), pntr(m+1)
      real*8 val(*), x(n), y(m), yj
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer i,j,jbgn, jend, indxi, rem, chunk, k, incn
      real*8 vali, yj1, yj2, yj3, yj4, yj5

*
*     --------------------------
*     First executable statement
*     --------------------------
*
      if (nrhs.eq.1) then
      do 110 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         yj1 = 0.0
         do 120 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incn = indx(i)
           yj1 = yj1 + val(i) * x(incn)
 120     continue
         incn = j
         y(incn) = yj1
 110    continue


      else if (nrhs.eq.2) then

      do 210 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         yj1 = 0.0
         yj2 = 0.0
         do 220 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incn = indx(i)
           yj1 = yj1 + val(i) * x(incn)
           incn = incn + n
           yj2 = yj2 + val(i) * x(incn)
 220     continue
         incn = j
         y(incn) = yj1
         incn = incn + n
         y(incn) = yj2
 210  continue

      else if (nrhs.eq.3) then

      do 310 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         yj1 = 0.0
         yj2 = 0.0
         yj3 = 0.0
         do 320 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incn = indx(i)
           yj1 = yj1 + val(i) * x(incn)
           incn = incn + n
           yj2 = yj2 + val(i) * x(incn)
           incn = incn + n
           yj3 = yj3 + val(i) * x(incn)
 320     continue
         incn = j
         y(incn) = yj1
         incn = incn + n
         y(incn) = yj2
         incn = incn + n
         y(incn) = yj3
 310  continue

      else if (nrhs.eq.4) then

      do 410 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         yj1 = 0.0
         yj2 = 0.0
         yj3 = 0.0
         yj4 = 0.0
         do 420 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incn = indx(i)
           yj1 = yj1 + val(i) * x(incn)
           incn = incn + n
           yj2 = yj2 + val(i) * x(incn)
           incn = incn + n
           yj3 = yj3 + val(i) * x(incn)
           incn = incn + n
           yj4 = yj4 + val(i) * x(incn)
 420     continue
         incn = j
         y(incn) = yj1
         incn = incn + n
         y(incn) = yj2
         incn = incn + n
         y(incn) = yj3
         incn = incn + n
         y(incn) = yj4
 410  continue

      else if (nrhs.eq.5) then

      do 510 j = 1, n
         jbgn = pntr(j)
         jend = pntr(j+1) - 1
         incn = 0
         yj1 = 0.0
         yj2 = 0.0
         yj3 = 0.0
         yj4 = 0.0
         yj5 = 0.0
         do 520 i = jbgn, jend
           indxi = indx(i)
           vali = val(i)
           incn = indx(i)
           yj1 = yj1 + val(i) * x(incn)
           incn = incn + n
           yj2 = yj2 + val(i) * x(incn)
           incn = incn + n
           yj3 = yj3 + val(i) * x(incn)
           incn = incn + n
           yj4 = yj4 + val(i) * x(incn)
           incn = incn + n
           yj5 = yj5 + val(i) * x(incn)
 520     continue
         incn = j
         y(incn) = yj1
         incn = incn + n
         y(incn) = yj2
         incn = incn + n
         y(incn) = yj3
         incn = incn + n
         y(incn) = yj4
         incn = incn + n
         y(incn) = yj5
 510  continue


      endif
      return
      end
