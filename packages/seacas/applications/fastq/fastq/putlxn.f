C $Id: putlxn.f,v 1.1 1990/11/30 11:13:59 gdsjaar Exp $
C $Log: putlxn.f,v $
C Revision 1.1  1990/11/30 11:13:59  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]PUTLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PUTLXN (MXND, NL, LXN, NUID, NODE, LINES, NAVAIL,
     &   IAVAIL, NNN, ERR, NOROOM)
C***********************************************************************
C
C  SUBROUTINE PUTLXN = DEFINE THE LINES FOR NODE AS  (LINES (I), I=1, NL)
C
C***********************************************************************
C
C  NOTE:
C     SAME CONTINUATION ENTRIES ARE USED AS ALREADY IN USE
C     FOR THIS NODE.
C     THIS NODE WILL BE FLAGGED AS A BOUNDARY NODE IF
C     LXN (2, NODE) .LT. 0   (NOT IF LINES (2) .LT. 0)
C
C***********************************************************************
C
      DIMENSION LINES (NL), LXN (4, MXND), NUID (MXND)
C
      LOGICAL ERR, BOUND, NOROOM
C
      BOUND = .FALSE.
C
      IF (LXN (2, NODE) .LT. 0)BOUND = .TRUE.
      NN = NODE
      NDONE = 0
C
C  FILL IN NEXT 3  (4 LAST TIME) NODES
C
  100 CONTINUE
      N4 = LXN (4, NN)
      NR = MIN0 (4, NL - NDONE)
      DO 110 I = 1, NR
         J = NDONE + I
         LXN (I, NN) = LINES (J)
  110 CONTINUE
C
C  CLEAR REMAINING PORTION
C
      IF (NR .LT. 4) THEN
         NZ = NR + 1
         DO 120 I = NZ, 4
            LXN (I, NN) = 0
  120    CONTINUE
      ENDIF
C
C  TAG BOUNDARY NODES
C
      IF (BOUND)LXN (2, NN) =  - LXN (2, NN)
C
C  TAG CONTINUATIONS
C
      IF (NDONE .GT. 1)LXN (1, NN) =  - LXN (1, NN)
      IF (NDONE + 4 .GE. NL) THEN
C
C  COLLECT GARBAGE
C
  130    CONTINUE
         IF (N4 .GE. 0) THEN
C
C  UPDATE NNN
C
  140       CONTINUE
            IF (LXN (1, NNN) .NE. 0) THEN
               RETURN
            ELSE
               NNN = NNN - 1
               GOTO 140
            ENDIF
         ENDIF
C
         NR =  - N4
         N4 = LXN (4, NR)
         LXN (1, NR) = 0
         LXN (2, NR) = 0
         LXN (3, NR) = 0
         LXN (4, NR) = IAVAIL
         IAVAIL = NR
         NAVAIL = NAVAIL + 1
         GOTO 130
      ENDIF
C
C  NEED ANOTHER LINE IN THE TABLE
C
      NDONE = NDONE + 3
      NEXTR = IABS (N4)
      IF (N4 .LT. 0) THEN
         LXN (4, NN) =  - NEXTR
         NN = NEXTR
         GOTO 100
      ENDIF
C
C  RESERVE A NEW LINE IN LXN TABLE
C
      IF (NAVAIL .LT. 1) THEN
         WRITE ( * , 10000)NODE
         ERR = .TRUE.
         NOROOM = .TRUE.
         RETURN
      ENDIF
      NEW = IAVAIL
      IF (NEW .GT. NNN)NNN = NEW
      IAVAIL = LXN (4, IAVAIL)
      NAVAIL = NAVAIL - 1
      LXN (4, NEW) = 0
      NUID (NEW) = 0
      NEXTR = NEW
C
C  INSERT LINK TO NEXT LINE
C
      LXN (4, NN) =  - NEXTR
      NN = NEXTR
      GOTO 100
C
10000 FORMAT (' NODE TABLE OVERFLOW IN PUTLXN - NODE = ', I5)
C
      END
