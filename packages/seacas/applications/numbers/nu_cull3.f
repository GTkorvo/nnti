C $Id: cull3.f,v 1.1 1991/02/21 15:42:47 gdsjaar Exp $
C $Log: cull3.f,v $
C Revision 1.1  1991/02/21 15:42:47  gdsjaar
C Initial revision
C
      SUBROUTINE CULL3 (DIRCOS, MASSLV, NIQM)
      DIMENSION DIRCOS(5,*), MASSLV(2,*)
C
      J = 0
      DO 10 I=1, NIQM
          IF (MASSLV(2,I) .NE. 0) THEN
              J = J + 1
              MASSLV(1, J) = MASSLV(1, I)
              MASSLV(2, J) = MASSLV(2, I)
              DIRCOS(1, J) = DIRCOS(1, I)
              DIRCOS(2, J) = DIRCOS(2, I)
              DIRCOS(3, J) = DIRCOS(3, I)
              DIRCOS(4, J) = DIRCOS(4, I)
              DIRCOS(5, J) = DIRCOS(5, I)
          END IF
   10 CONTINUE
C
      NIQM = J
      RETURN
      END
