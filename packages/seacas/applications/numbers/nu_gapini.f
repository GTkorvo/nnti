C $Id: gapini.f,v 1.3 1999/02/16 21:38:00 gdsjaar Exp $
C $Log: gapini.f,v $
C Revision 1.3  1999/02/16 21:38:00  gdsjaar
C Converted to read exodusII database format.  Somewhat tested, not
C ready for production yet.
C
C Revision 1.2  1997/06/20 19:11:26  caforsy
C Port to ibm
C
C Revision 1.1.1.1  1991/02/21 15:43:14  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:43:13  gdsjaar
c Initial revision
c
      SUBROUTINE GAPINI (A, COORD, IDESS, NEESS, NNESS, IPEESS, IPNESS,
     *   LTEESS, LTNESS, FACESS, DISP, NUMNP, NDIM, NUMESS,
     *   TIME, ITMSEL, TITLE, IMAS, ISLV, DMAX, GMTHD)
C
      DIMENSION A(*), COORD(NUMNP,*), IDESS(*), NEESS(*),
     *   NNESS(*), IPEESS(*), IPNESS(*), LTEESS(*), LTNESS(*),
     *   FACESS(*), TIME(*), DISP(NUMNP,*)
      LOGICAL ITMSEL(*)
      CHARACTER*80 TITLE, STRA
      CHARACTER*8  GMTHD
      LOGICAL ERROR

      IFLGM = LOCINT (IMAS, NUMESS, IDESS)
      IFLGS = LOCINT (ISLV, NUMESS, IDESS)
C
      ERROR = .FALSE.
      IF (IFLGM .EQ. 0) THEN
         WRITE (STRA, 10) 'Master', IMAS
         CALL SQZSTR (STRA, LSTR)
         CALL PRTERR ('ERROR', STRA(:LSTR))
   10    FORMAT (1X,A,' Surface Flag ',I5,' not found. ')
         ERROR = .TRUE.
      END IF
      IF (IFLGS .EQ. 0) THEN
         WRITE (STRA, 10) 'Slave', ISLV
         CALL SQZSTR (STRA, LSTR)
         CALL PRTERR ('ERROR', STRA(:LSTR))
         ERROR = .TRUE.
      END IF
      IF (ERROR) RETURN
C
      NSEGM = NEESS(IFLGM)
      IPTRM = IPNESS(IFLGM)
C
      NSEGS = NEESS(IFLGS)
      IPTRS = IPNESS(IFLGS)
C
      MULT = 2 * NDIM - 2
      CALL MDRSRV ('MAPMAS', IMPMS, MULT*NSEGM)
      CALL MDRSRV ('MAPSLV', IMPSL, MULT*NSEGS)
      CALL MDRSRV ('ITEMP',  ITMP,  MAX(NUMNP,3*NSEGM))
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) THEN
         CALL MEMERR
         STOP
      END IF
C
      CALL UNIQUE (LTNESS(IPTRM), MULT*NSEGM, A(IMPMS), A(ITMP),
     *   NIQM, NUMNP)
      CALL MDRSRV ('MASSLV', IMSLV, 2*NIQM)
      CALL MDRSRV ('DIRCOS', IDCOS, (NDIM+2)*NIQM)
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) THEN
         CALL MEMERR
         STOP
      END IF
      CALL TRANIQ (LTNESS(IPTRM), A(IMPMS), A(IMSLV), MULT*NSEGM, 2)
C
      CALL UNIQUE (LTNESS(IPTRS), MULT*NSEGS, A(IMPSL), A(ITMP),
     *   NIQS, NUMNP)
      CALL MDRSRV ('NIQSLV', INQS, NIQS)
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) THEN
         CALL MEMERR
         STOP
      END IF
      CALL TRANIQ (LTNESS(IPTRS), A(IMPSL), A(INQS), MULT*NSEGS, 1)
C
      DMAX = DMAX**2
      IF (DMAX .EQ. 0.0) DMAX = 1.0E38

      IF (MAX(NUMNP, 3*NSEGM) .LT. 4*NIQS) THEN
         CALL MDLONG ('ITEMP', ITMP, 4*NIQS)
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP
         END IF
      END IF
      IF (NDIM .EQ. 3) THEN
         CALL DONRM3 (COORD, LTNESS(IPTRM), A(IMPMS), A(IDCOS),
     *      A(ITMP), NSEGM, NIQM, NUMNP)
         IF (GMTHD .EQ. 'DISTANCE') THEN
            CALL GMDIS3 (COORD, A(IDCOS), A(IMSLV), A(INQS), A(ITMP),
     *         A(ITMP+3*NIQS), NIQM, NIQS, DMAX, NUMNP)
         ELSE IF (GMTHD .EQ. 'NORMAL') THEN
            CALL GMTCH3 (COORD, A(IDCOS), A(IMSLV), A(INQS), A(ITMP),
     *         A(ITMP+3*NIQS), NIQM, NIQS, DMAX, NUMNP)
         END IF
      ELSE
         CALL DONRM2 (COORD, LTNESS(IPTRM), A(IMPMS), A(IDCOS),
     *      A(ITMP), NSEGM, NIQM, NUMNP)
         IF (GMTHD .EQ. 'DISTANCE') THEN
            CALL GMDIS2 (COORD, A(IDCOS), A(IMSLV), A(INQS), A(ITMP),
     *         NIQM, NIQS, DMAX, NUMNP)
         ELSE IF (GMTHD .EQ. 'NORMAL') THEN
            CALL GMTCH2 (COORD, A(IDCOS), A(IMSLV), A(INQS), A(ITMP),
     *         NIQM, NIQS, DMAX, NUMNP)
         END IF
      END IF
      CALL GAPOUT (A(IDCOS), A(IMSLV), NIQM, NDIM, IDESS(IFLGM),
     *   IDESS(IFLGS), GMTHD)

      CALL MDDEL ('MAPMAS')
      CALL MDDEL ('MAPSLV')
      CALL MDDEL ('ITEMP' )
      CALL MDDEL ('MASSLV')
      CALL MDDEL ('DIRCOS')
      CALL MDDEL ('NIQSLV')
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) THEN
         CALL MEMERR
         STOP
      END IF
      RETURN
      END
