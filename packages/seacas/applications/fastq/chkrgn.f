C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C $Id: chkrgn.f,v 1.4 2004/01/21 05:18:39 gdsjaar Exp $
C $Log: chkrgn.f,v $
C Revision 1.4  2004/01/21 05:18:39  gdsjaar
C Initialized several variables identified by valgrind.
C
C Revision 1.3  2001/11/05 13:26:51  gdsjaar
C  Fixed array boundary problem in region check code.
C
C Revision 1.2  1990/11/30 11:25:08  gdsjaar
C Added initialization and MDSTAT calls
C
c Revision 1.1.1.1  90/11/30  11:04:38  gdsjaar
c FASTQ Version 2.0X
c 
c Revision 1.1  90/11/30  11:04:37  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]CHKRGN
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO CHKHOL TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE CHKRGN (IA, L, MP, ML, MS, MR, MSC, N24, IPOINT, COOR,
     &   IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN,
     &   ISIDE, NLPS, IFLINE, ILLIST, IREGN, NSPR, IFSIDE, ISLIST,
     &   NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB,
     &   IFHOLE, NHPR, IHLIST, LINKP, LINKL, LINKS, LINKR, LINKSC,
     &   LINKPB, LINKLB, LINKSB, RSIZE, SCHEME, DEFSCH, NPREGN, NPSBC,
     &   NPNODE, MAXNP, MAXNL, MAX3, X, Y, NID, LISTL, NNPS, ANGLE,
     &   MARKED, MXND, MXNPER, MXNL, MAXNBC, MAXSBC, AMESUR, XNOLD,
     &   YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &   NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &   IDIVIS, SIZMIN, EMAX, EMIN, NOROOM, ERRCHK, ERR)
C***********************************************************************
C
C  SUBROUTINE CHKRGN - CHECK THAT A REGION MAY BE MESHED
C
C***********************************************************************
C
      DIMENSION IA(1)
      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML), ISIDE(MS), NLPS(MS)
      DIMENSION IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION SCHEME(MSC), RSIZE (MR)
      DIMENSION NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKR(2, MR), LINKSC(2, MR), LINKPB(2, MP)
      DIMENSION LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION X(MAXNP), Y(MAXNP), NID(MAXNP)
      DIMENSION LISTL(MAXNL), MARKED(3, MAXNL)
      DIMENSION NNPS(MAX3), ANGLE(MAXNP)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
C
      DIMENSION IDUMMY(1)
C
      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      CHARACTER*72 SCHEME, DEFSCH, SCHSTR
C
      LOGICAL NOROOM, EVEN, ERR, NORM, CCW, REAL, ADDLNK, REMESH
      LOGICAL PENTAG, TRIANG, TRNSIT, HALFC, COUNT, ERRCHK
C
      ipntr = 0
      addlnk = .false.
      COUNT = .TRUE.
      IF (REMESH) THEN
         EVEN = .TRUE.
      ELSE
         EVEN = .FALSE.
      ENDIF
      REAL = .FALSE.
C
C  CHECK TO MAKE SURE CONNECTING DATA FOR THE REGION EXISTS
C  AND FILL IN ANY BLANK INTERVALS ACCORDING TO THE GIVEN SIZE
C  FOR THE REGION AND THE LINE'S LENGTH
C
      CALL DATAOK (MP, ML, MS, MR, L, IREGN(L), COOR, ILINE, LTYPE,
     &   NINT, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, RSIZE(L), ERRCHK, ERR)
      IF (ERR) THEN
         WRITE (*, 10000) IREGN(L)
         ADDLNK = .TRUE.
         IMINUS = -L
         CALL LTSORT (MR, LINKR, IREGN(L), IMINUS, ADDLNK)
         ADDLNK = .FALSE.
C
C  CALCULATE THE PERIMETER OF THE REGION
C
      ELSE
         KNBC = 0
         KSBC = 0
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO PERIM TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
         CALL PERIM (MP, ML, MS, NSPR(L), MAXNL, MAXNP, 1, 1, KNBC,
     &      KSBC, IREGN(L), IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT,
     &      FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST,
     &      ISLIST(IFSIDE(L)), NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB,
     &      NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS, LINKPB, LINKLB,
     &      LINKSB, X, Y, NID, NPER, LISTL, NL, IDUMMY, MARKED, EVEN,
     &      REAL, ERR, CCW, COUNT, NOROOM, AMESUR, XNOLD, YNOLD, NXKOLD,
     &      MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD,
     &      NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &      IDIVIS, SIZMIN, EMAX, EMIN)
         IF ((NPER .LE. 0) .OR. (ERR)) THEN
            WRITE (*, 10010) IREGN(L)
            ADDLNK = .TRUE.
            IMINUS = -L
            CALL LTSORT (MR, LINKR, IREGN(L), IMINUS, ADDLNK)
            ADDLNK = .FALSE.
            GO TO 120
         END IF
C
C  WHEN CHECKING THE MAXIMUMS - ADD ENOUGH FOR ONE MORE INTERVAL
C  ON THE LINE AS THIS LINE MAY BE INCREMENTED BY ONE IF THE
C  PERIMETER IS ODD
C
         MAXNBC = MAX(MAXNBC, KNBC + 3)
         MAXSBC = MAX(MAXSBC, KSBC + 3)
         MXNL = MAX(MXNL, NL)
C
C  GET THE REGION SCHEME
C
         CALL LTSORT (MR, LINKSC, ABS(IREGN(L)), IPNTR, ADDLNK)
         IF ((IREGN(L) .LE. N24) .AND. (IPNTR .GT. 0)) THEN
            SCHSTR = SCHEME(IPNTR)
         ELSE
            SCHSTR = DEFSCH
         END IF
C
C  SEE IF A TRIANGULAR, PENTAGON, SEMICIRCLE, OR A TRANSITION
C  REGION HAS BEEN FLAGGED
C
         PENTAG = .FALSE.
         TRNSIT = .FALSE.
         TRIANG = .FALSE.
         CALL STRCUT (SCHSTR)
         CALL STRLNG (SCHSTR, LENSCH)
         DO 100 J = 1, LENSCH
            IF ((SCHSTR(J:J) .EQ. 'T') .OR.
     &         (SCHSTR(J:J) .EQ. 't')) THEN
               IF (NPER .GE. 6) THEN
                  TRIANG = .TRUE.
               ELSE
                  CALL MESAGE ('TRIANGULAR REGION MESH NOT')
                  CALL MESAGE ('POSSIBLE WITH PERIMETER < 6')
                  CALL MESAGE ('REGULAR PROCESSING ASSUMED')
               END IF
               GO TO 110
            ELSE IF ((SCHSTR(J:J) .EQ. 'U') .OR.
     &         (SCHSTR(J:J) .EQ. 'u')) THEN
               IF (NPER .GE. 10) THEN
                  PENTAG = .TRUE.
               ELSE
                  CALL MESAGE ('PENTAGON REGION MESH NOT')
                  CALL MESAGE ('POSSIBLE WITH PERIMETER < 10')
                  CALL MESAGE ('REGULAR PROCESSING ASSUMED')
               END IF
               GO TO 110
            ELSE IF ((SCHSTR(J:J) .EQ. 'B') .OR.
     &         (SCHSTR(J:J) .EQ. 'b')) THEN
               IF (NPER .GE. 8) THEN
                  TRNSIT = .TRUE.
                  HALFC = .FALSE.
               ELSE
                  CALL MESAGE ('TRANSITION REGION GENERATION NOT')
                  CALL MESAGE ('POSSIBLE WITH NO. IN PERIMETER < 8')
                  CALL MESAGE ('REGULAR PROCESSING WILL BE ATTEMPTED')
               END IF
               GO TO 110
            ELSE IF ((SCHSTR(J:J) .EQ. 'C') .OR.
     &         (SCHSTR(J:J) .EQ. 'c'))
     &         THEN
               IF (NPER .GE. 8) THEN
                  TRNSIT = .TRUE.
                  HALFC = .TRUE.
               ELSE
                  CALL MESAGE ('SEMICIRCLE REGION GENERATION NOT')
                  CALL MESAGE ('POSSIBLE WITH NO. IN PERIMETER < 8')
                  CALL MESAGE ('REGULAR PROCESSING WILL BE ATTEMPTED')
               END IF
               GO TO 110
            END IF
  100    CONTINUE
  110    CONTINUE
C
C  SET UP THE TRIANGLE DIVISIONS, AND FIND THE CENTER POINT
C
         IF (TRIANG) THEN
            CALL GETM3 (ML, MS, MAX3, NSPR(L), ISLIST(IFSIDE(L)), NINT,
     &         IFLINE, NLPS, ILLIST, LINKL, LINKS, X, Y, NID, NNPS,
     &         ANGLE, NPER, M1A, M1B, M2A, M2B, M3A, M3B, XCEN, YCEN,
     &         CCW, ERR)
C
C  CHECK FOR MAXIMUM DIMENSIONS NEEDED FOR EACH REGION
C  ASSUMING THAT 10 NECKLACES WILL BE ADEQUATE
C
            MXTEST = ((M1A + 1)*(M3B + 1)) + ((M1B + 1)*(M2A + 1))
     &         + ((M2B + 1)*(M3A + 1)) + (10*(NPER + 1)) + (NPER*2)
            MXND = MAX(MXTEST, MXND)
            MXNPER = MAX(MXNPER, (NPER + 2) * 2)
C
C  SET UP THE TRANSITION DIVISIONS, AND FIND THE CENTER POINT
C
         ELSE IF (TRNSIT) THEN
            CALL GETTRN (ML, MS, MAX3, NSPR(L), ISLIST(IFSIDE(L)), NINT,
     &         IFLINE, NLPS, ILLIST, LINKL, LINKS, X, Y, NID, NNPS,
     &         ANGLE, NPER, I1, I2, I3, I4, I5, I6, I7, I8, XCEN1,
     &         YCEN1, XCEN2, YCEN2, XMID1, YMID1, XMID2, YMID2, CCW,
     &         HALFC, ERR)
C
C  CHECK FOR MAXIMUM DIMENSIONS NEEDED FOR EACH REGION
C  ASSUMING THAT 10 NECKLACES WILL BE ADEQUATE
C
            MXTEST = ((I2 - I1)*(NPER - I8))
     &         + ((I3 - I2)*(NPER - I8)) + ((I3 - I2)*(I2 - I2))
     &         + ((I4 - I3)*(I7 - I6)) + ((I5 - I4)*(I6 - I5))
     &         + ((I5 - I4)*(I7 - I6)) + (10*(NPER + 1)) + (NPER*2)
            MXND = MAX(MXTEST, MXND)
            MXNPER = MAX(MXNPER, (NPER + 2) * 2)
C
C  SET UP THE PENTAGON DIVISIONS, AND FIND THE CENTER POINT
C
         ELSE IF (PENTAG) THEN
            CALL GETM5 (IA, ML, MS, MAX3, NSPR(L), ISLIST(IFSIDE(L)),
     &         NINT, IFLINE, NLPS, ILLIST, LINKL, LINKS, X, Y, NID,
     &         NNPS, ANGLE, NPER, M1A, M1B, M2, M3A, M3B, M4A, M4B,
     &         M5, MC, XCEN, YCEN, CCW, ERR)
C
C  CHECK FOR MAXIMUM DIMENSIONS NEEDED FOR THE REGION
C  ASSUMING THAT 10 NECKLACES WILL BE ADEQUATE
C
            MXTEST = (M1B*M2) + (M4A*M3B) + (M4B*M5) + (10*(NPER + 1))
            MXND = MAX(MXTEST, MXND)
            MXNPER = MAX(MXNPER, (NPER + 2) *2)
C
C  CALCULATE THE BASE OF THE RECTANGLE FOR THE REGION
C
         ELSE
            CALL GETM1 (ML, MS, MAX3, NSPR(L), ISLIST(IFSIDE(L)), NINT,
     &         IFLINE, NLPS, ILLIST, LINKL, LINKS, X, Y, NID, NNPS,
     &         ANGLE, NPER, SCHSTR, M1, CCW, NORM, REAL, ERR)
            IF (ERR) THEN
               WRITE (*, 10020) IREGN(L)
               ADDLNK = .TRUE.
               IMINUS = -L
               CALL LTSORT (MR, LINKR, IREGN(L), IMINUS, ADDLNK)
               ADDLNK = .FALSE.
               GO TO 120
            END IF
            M2 = NPER/2 - M1
C
C  CHECK FOR MAXIMUM DIMENSIONS NEEDED FOR EACH REGION
C  ASSUMING THAT 10 NECKLACES WILL BE ADEQUATE
C
            MXTEST = ((M1 + 1)*(M2 + 1)) + (10*(M1 + M2 + 2))
            MXND = MAX(MXTEST, MXND)
            MXNPER = MAX(MXNPER, NPER + 4)
         END IF
C
C  FLAG THE REGION AS BEING PROCESSABLE
C
         IREGN(L) = -IREGN(L)
C
C  MARK THE LINES AND POINTS IN THE REGION AS BEING USED
C
         CALL MKUSED (MAXNL, MP, ML, LISTL, IPOINT, NINT, LINKP, LINKL,
     &      LCON, NL)
C
C  CHECK ALL THE HOLES IN THE REGION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO CHKHOL TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
         CALL CHKHOL (IA, L, MP, ML, MS, MR, MSC, IPOINT, COOR,
     &      IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN,
     &      ISIDE, NLPS, IFLINE, ILLIST, IREGN, NSPR, IFSIDE, ISLIST,
     &      NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB,
     &      IFHOLE, NHPR, IHLIST, LINKP, LINKL, LINKS, LINKR, LINKSC,
     &      LINKPB, LINKLB, LINKSB, RSIZE, NPREGN, NPSBC, NPNODE,
     &      MAXNP, MAXNL, MXNPER, KNBC, KSBC, X, Y, NID, LISTL, MARKED,
     &      MXNL, MAXNBC, MAXSBC, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD,
     &      LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD, NNXK,
     &      REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN,
     &      EMAX, EMIN, NOROOM, ERRCHK, ERR)

      END IF
C
  120 CONTINUE
      RETURN
C
10000 FORMAT (' ** ERROR - DATA PROBLEMS FOR REGION:', I5, ' **')
10010 FORMAT (' ** ERROR - PERIMETER GENERATION ERRORS FOR REGION:'
     &   , I5, ' **')
10020 FORMAT (' ** ERROR - MAPPING BASE GENERATION ERRORS FOR REGION:'
     &   , I5, ' **')
      END
