C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
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

C $Id: pltpal.f,v 1.1 1993/07/16 16:49:05 gdsjaar Exp $ 
C $Log: pltpal.f,v $
C Revision 1.1  1993/07/16 16:49:05  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTPAL(COL,R,G,B)
      CHARACTER*10 ECOLOR
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)

      IF (COL.LT.8. .OR. COL.GT.15.) THEN
         CALL CHRIC(INT(COL),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Illegal palette color '//ECOLOR(1:L)//
     *               ' passed to PLTPAL; range is 8-15.',2)
         RETURN

      END IF

      IF (R.LT.0. .OR. R.GT.1.) THEN
         CALL CHRIC(INT(R),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Red value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
      END IF

      IF (G.LT.0. .OR. G.GT.1.) THEN
         CALL CHRIC(INT(G),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Green value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
      END IF

      IF (B.LT.0. .OR. B.GT.1.) THEN
         CALL CHRIC(INT(B),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Blue value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
      END IF

      PALETT(1,INT(COL)) = R
      PALETT(2,INT(COL)) = G
      PALETT(3,INT(COL)) = B
      CALL VDSTCO(1,INT(COL),PALETT(1,INT(COL)),0)
      RETURN

      END
