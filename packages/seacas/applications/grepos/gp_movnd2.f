C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
      subroutine movnd2(numnp, ndim, x, y, vnorm, 
     *  neessm, nnessm, ltnesm, toler, delmax, index, vector, gap)
C=======================================================================
      
      REAL X(*), Y(*), VNORM(3,*)
      DIMENSION LTNESM(2,*)
      integer index(*)
      REAL VECTOR(2)
      
      dmax = 0.0
      dmin = 1.0e15
      match = 0

      do 130 inod = 1, numnp
        pmin = 1.0e38
        if (vnorm(1,inod) .ne. 0.0 .or. vnorm(2,inod) .ne. 0.0) then
          
          X0 = X(inod)
          Y0 = Y(inod)
          
          AI = VNORM(1, inod)
          BJ = VNORM(2, inod)
C ... Node movement (delta) = (xnew-X0)**2 + (ynew-Y0)**2
C                           = (x0+t*ai-x0)**2 + (y0+t*bj-y0)**2
C                           = t**2 * (ai**2 + bj**2)
C    Want delta < delmax  ==> t**2 * (ai**2 + bj**2) < delmax**2
C                         ==> t**2 < delmax**2 / (ai**2 + bj**2) = tmax
          
          tmax = delmax**2 / (ai**2 + bj**2)
          
          do 110 iseg = 1, neessm
            XI = x(LTNESM(1,ISEG))
            YI = y(LTNESM(1,ISEG))

            XJ = x(LTNESM(2,ISEG))
            YJ = y(LTNESM(2,ISEG))

C ... If denom == 0, then node normal is parallel to plane
            denom = (yj-yi)*ai - (xj-xi)*bj
            if (denom .ne. 0.0) then
              T = ((xj-xi)*(y0-yi) - (yj-yi)*(x0-xi))/denom
              S = (    ai *(y0-yi) -     bj *(x0-xi))/denom

              if (t .ge. 0.0 .and. t**2 .le. tmax .and.
     *          0.0-toler .le. S .and. S .le. 1.0+toler) then
C ... If we made it this far, then the intersection point is inside the
C     face. Save the minimum distance found so far.
                delta = t**2 * (ai**2 + bj**2)
                dmin = min(delta, dmin)
                dmax = max(delta, dmax)
                match = match + 1

                go to 120
              else
C ... The node is outside the tolerance, 
C     for this face/node combination. Save the minimum for all face/node comb
                if (S .lt. 0.0) then
                  S = -S
                else if (S .gt. 1.0) then
                  S = S - 1.0
                endif
                
                pmin = min(S, pmin)
              end if
            end if
 110      continue
        end if
 120    continue
 130  continue
      
C ... Update the node positions based on the minimum distance found
C     and the specified vector.
      if (match .gt. 0) then
        dmin = dmin - gap
        AI = vector(1)
        BJ = vector(2)

        do 140 inod=1, numnp
          if (index(inod) .eq. 1) then
            X0 = X(inod)
            Y0 = Y(inod)

C ... Update the nodes position (Currently, assumes in vector direction)
            X(inod) = X0 + dmin * AI
            Y(inod) = Y0 + dmin * BJ
          end if
 140    continue
        write (*, 10020) dmin
      else
        write (*,*) 'No node movement.'
      end if
      
10020 format(/,'Node movement = ',1pe11.4)
      return
      end
      
