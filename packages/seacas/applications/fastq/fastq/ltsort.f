C $Id: ltsort.f,v 1.3 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: ltsort.f,v $
C Revision 1.3  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1999/06/17 19:02:22  gdsjaar
C Fixed several problems related to holes.  In several places, a
C nonpositive integer was being used to index into an array.  This seems
C to fix all of those cases.  I'm not sure if I fixed the true cause of
C these errors or just the symptom though...
C
C Revision 1.1.1.1  1990/11/30 11:11:44  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:11:42  gdsjaar
c Initial revision
c 
C

CC* FILE: [.MAIN]LTSORT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LTSORT (MDIM, LINK, ID, IPNTR, ADDLNK)
C***********************************************************************
C
C  SUBROUTINE LTSORT = LOOKUP TABLE SORT FOR DATA POINTER ARRAYS
C
C***********************************************************************
C
C  VARIABLES USED:
C     MDIM   = DIMENSION OF LINK ARRAY,  AND BASE FOR LOOKUP START
C     LINK   = LOOKUP TABLE ARRAY OF ID'S AND POINTERS
C              LINK (1, I) = ID VALUE STORED IN I'TH ROW  (0 IF EMPTY)
C              LINK (2, I) = DATA POINTER ASSOCIATED W/THIS I'TH ID VALUE
C     ID     = THE ID OF THE DATA BEING FOUND OR PLACED
C     IPNTR  = THE DATA POINTER ASSOCIATED WITH THE ID BEING USED
C     ADDLNK = .TRUE. IF DATA IS BEING PLACED IN THE LOOKUP TABLE
C            = .FALSE. IF DATA IS BEING FOUND ONLY
C
C***********************************************************************
C
      DIMENSION LINK (2, MDIM)
C
      LOGICAL ADDLNK
C
C  CALCULATE THE BEGINNING LOOKUP VALUE
C
      if (id .lt. 0) stop 'LTSORT: Internal error'
      
      HOLD = FLOAT (ID) * 3.1830989
      LOOKUP =  INT((HOLD - IFIX (HOLD)) * FLOAT (MDIM) + 1)
C
C  SET UP THE LOOP TO ONLY SEARCH THROUGH THE TABLE ONCE
C
      DO 100 I = 1, MDIM
C
C  IF LOOKUP SPOT IS EMPTY THEN FILL AND RETURN IF ADDING AND IPNTR .NE. 0
C  OR FLAG IPNTR AS BEING EMPTY AND RETURN IF FINDING
C
         IF (LINK (1, LOOKUP) .EQ. 0) THEN
            IF ( (ADDLNK) .AND. (IPNTR .NE. 0)) THEN
               LINK (1, LOOKUP) = ID
               LINK (2, LOOKUP) = IPNTR
            ELSEIF (IPNTR .NE. 0) THEN
               IPNTR = 0
            ENDIF
            RETURN
C
C  IF LOOKUP SLOT IS FULL,  CHECK TO SEE IF IT MATCHES THE CURRENT ID
C  IF IT MATCHES AND IF ADDING,  SET THE NEW POINTER  (OVERWRITE)
C  IF IT MATCHES AND IF FINDING,  RETURN THE CORRECT POINTER
C  IF NO MATCH,  THEN INCREMENT LOOKUP AND TRY AGAIN IN THE TABLE
C
         ELSE
            IF (ID .EQ. LINK (1, LOOKUP)) THEN
               IF (ADDLNK) THEN
                  LINK (2, LOOKUP) = IPNTR
               ELSE
                  IPNTR = LINK (2, LOOKUP)
               ENDIF
               RETURN
            ELSE
               LOOKUP = LOOKUP + 1
               IF (LOOKUP .GT. MDIM)LOOKUP = 1
            ENDIF
         ENDIF
  100 CONTINUE
C
C  ACT ON THE EXHAUSTED SEARCH
C
      IF (ADDLNK) THEN
         CALL MESAGE ('LOOKUP TABLE OVERFLOW')
         CALL MESAGE ('SERIOUS DATA PROBLEMS HAVE BEEN CAUSED')
      ELSE
         IPNTR = 0
      ENDIF
      RETURN
C
      END
