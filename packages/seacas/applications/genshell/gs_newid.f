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

      SUBROUTINE NEWID (TYPE, IDLST, NUMID, IDNEW, IDOLD)
C
      CHARACTER*(*) TYPE
      DIMENSION IDLST(*)
      CHARACTER*80 STRING
      CHARACTER*8 STRA

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in NEWID')
         RETURN
      END IF

C ... Check for existance of new id in old list - not allowed

      IF (NUMID .LE. 0) RETURN

      IMAT = LOCINT (IDNEW, NUMID, IDLST)

      IF (IMAT .NE. 0) THEN
         CALL PRTERR ('ERROR', 'Cannot change to an existing ID')
         RETURN
      END IF

C ... Determine location of ID to be changed

      IMAT = LOCINT (IDOLD, NUMID, IDLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) STRA, IDOLD
   90    FORMAT (A,1X,I5,' does not exist')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         RETURN
      ELSE
         IDLST(IMAT) = IDNEW
         WRITE (STRING, 100) STRA, IDOLD, STRA, IDNEW
  100    FORMAT (A,1X,I6,' changed to ',A,1X,I6)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
      END IF

      RETURN
      END
