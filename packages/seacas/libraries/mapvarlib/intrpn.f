C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
C retains certain rights in this software.
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

C=======================================================================
*DECK,INTRPN
      SUBROUTINE INTRPN(ICONA,SOLNA,IELPT,STRPT,
     &                  SOLNB,NDLSTB,XB,YB,ZB,
     &                  IDBLK,TIMES,INSUB,SN)
C     
C     ******************************************************************
C     
C     SUBROUTINE TO CONTROL INTERPOLATION OF NODAL RESULTS FROM MESH-A
C     TO MESH-B
C     INTERPOLATED SOLUTION IS WRITTEN TO MESH-C EXODUS FILE
C
C     Calls subroutine SHAPEF, ININOD
C
C     Called by MAPVAR
C     
C     ******************************************************************
C
C ICONA   INT   Connectivity of donor Mesh (1:nelnda,1:numeba)
C SOLNA   REAL  Nodal variables  for donor mesh
C IELPT   INT   The element in donor mesh within which the point 
C               (node in recipient mesh) is found
C STRPT   REAL  The isoparametric coords of point in IELPT element
C SOLNB   REAL  Nodal variables for recipient mesh
C NDLSTB  INT   List of recipient mesh nodes in current element block
C SOLN    REAL  SOLNA vector local to each donor mesh element
C TIMES   REAL  Array of times - just passed through
C INSUB   INT   Number of times into subroutine for this element block
C               1-first time; >1-second,etc time;
C               used to control mapping for many element blocks to one
C
C     ******************************************************************
C
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'aexds1.blk'
      include 'contrl.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'steps.blk'
      include 'tapes.blk'
      include 'varnpt.blk'
C
      DIMENSION XB(*), YB(*), ZB(*), TIMES(*)     
      DIMENSION ICONA(NELNDA,*), SOLNA(NODESA,NVARNP)
      DIMENSION SOLNB(NODESB,NVARNP), NDLSTB(*)
      DIMENSION IELPT(*), STRPT(3,NODESB)
      DIMENSION SOLN(27), SN(*)
C     
C     ******************************************************************
C
C set up time steps
C
      IF (ISTEP .EQ. -1)THEN
        NTM = NTIMES
      ELSE
        NTM = 1
      END IF
C
      DO 5 IST = 1, NTM
        IF (ISTEP .EQ. -1)THEN
          ISTP = IST
        ELSE
          ISTP = ISTEP
        END IF
C     
C Start interpolation
C
        DO 10 IVAR = 1, NVARNP
C
C For IDEF = 2 do mesh annealing (used by GOMA); write out all
C displacements as zero.
C
          IF (IDEF .EQ. 2 .AND. (IVAR .EQ. IXDIS .OR. IVAR .EQ. IYDIS
     &        .OR. IVAR .EQ. IZDIS))GO TO 10
C
C
C If first time into INTRPN for this element block, initialize
C else you are mapping many to one and retrieve partially mapped
C results from temporary storage in EXODUS
C
          IF (INSUB .EQ. 1)THEN
            CALL ININOD(SOLNB,IVAR,TIMES,ISTP,IDBLK,NDLSTB,XB,YB,ZB,
     &                  SN)
          ELSE
            CALL EXGNV(NTP4EX,IST,IVAR,NODESB,SOLNB(1,IVAR),IERR)
          END IF
C
C Get nodal results on donor mesh
C
          CALL EXGNV(NTP2EX,ISTP,IVAR,NODESA,SOLNA(1,IVAR),IERR)
C     
C Loop on nodes in recipient mesh
C     
          DO 30 I = 1,NUMNDB
            NEL = IELPT(I)
            IF (NEL .NE. 0) THEN
C     
C Set parameters for element in donor mesh
C     
              S = STRPT(1,I)
              T = STRPT(2,I)
              R = 0.
              IF (NDIMB.EQ.3) R = STRPT(3,I)
C              NNODES = NNELM(ITYPE)
              NNODES = NELNDA
              IF (ITYPE .EQ. 6) NNODES = 4
              DO 20 J = 1,NNODES
                INODE = ICONA(J,NEL)
                SOLN(J) = SOLNA(INODE,IVAR)
 20           CONTINUE
C     
C Shape function
C
              CALL SHAPEF(ITYPE,S,T,R,SOLN,BVALUE)
              SOLNB(NDLSTB(I),IVAR) = BVALUE
            END IF
 30       CONTINUE
C
C Save results, it doesn't matter if they are preliminary or final
C
          CALL EXPNV(NTP4EX,IST,IVAR,NODESB,SOLNB(1,IVAR),IERR)

 10     CONTINUE
  5   CONTINUE
C
      RETURN
      END
