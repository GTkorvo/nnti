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

C $Id: lxrem.f,v 1.1 1993/07/16 16:46:46 gdsjaar Exp $ 
C $Log: lxrem.f,v $
C Revision 1.1  1993/07/16 16:46:46  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE LXREM(LINE,L)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER*80 TMPLIN

      K = INDEX(ILINE(JLINE:),CHAR(0))
      IF (K.EQ.0) THEN
         JLINE = 504
         L = 0
         LINE = ' '

      ELSE
         L = K - 1
         IF (L.GT.LEN(LINE)) THEN
            L = LEN(LINE)
            LINE = ILINE(JLINE:JLINE+L-1)
            TMPLIN = 'Remainder truncated:'//LINE
            CALL LXERR(TMPLIN,1)

         ELSE IF (L.GT.0) THEN
            LINE = ILINE(JLINE:JLINE+K-1)
         END IF

         JLINE = JLINE + K
      END IF

      RETURN

      END
