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
      subroutine centroid(x, y , z, xcent, ycent, zcent,
     *  nelblk, numelb, numlnk, link, ndim)
      real x(*), y(*), z(*)
      real xcent(*), ycent(*), zcent(*)
      integer numelb(*), numlnk(*), link(*)

      IELNK = 0
      IE = 0
      DO 100 IELB = 1, nelblk
         IS = IE + 1
         IE = IE + NUMELB(IELB)
         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)

         CALL cent1(x, y, z, xcent(is), ycent(is), zcent(is),
     *     NUMELB(IELB), NUMLNK(IELB), LINK(ISLNK), NDIM)
  100 CONTINUE

      RETURN
      END
      

      subroutine cent1(x, y, z, xcent, ycent, zcent,
     *  numelb, numlnk, link, ndim)

      real x(*), y(*), z(*)
      real xcent(*), ycent(*), zcent(*)
      integer numelb, numlnk
      integer link(numlnk,*)
      integer ndim

      do 20 ne=1, numelb
        xc = 0.0
        yc = 0.0
        zc = 0.0
        do 10 j=1, numlnk 
          xc = xc + x(link(j,ne))
          yc = yc + y(link(j,ne))
          if (ndim .eq. 3) then
            zc = zc + z(link(j,ne))
          end if
 10     continue
        
        rnodes = numlnk
        xcent(ne) = xc / rnodes
        ycent(ne) = yc / rnodes
        if (ndim .eq. 3) then
          zcent(ne) = zc / rnodes
        end if
        
 20   continue
      return
      end
