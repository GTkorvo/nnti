C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C=======================================================================
      SUBROUTINE PRINIT (OPTION, NOUT, NDB, FILNAM, TITLE, NDIM, NUMNP,
     &     NUMEL, NELBLK, NUMNPS, LNPSNL, LNPSDF, NUMESS, 
     &     LESSEL, LESSNL, NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
C=======================================================================

C   --*** PRINIT *** (BLOT) Display database initial variables
C   --   Written by Amy Gilkey - revised 01/14/88
C   --
C   --PRINIT displays the database initial variables.
C   --
C   --Parameters:
C   --   OPTION  - IN - '*' to print all, else print options:
C   --                  'N' to print database name
C   --                  'T' to print title
C   --                  'I' to print number of nodes, etc.
C   --                  'S' to print node and side set information
C   --                  'V' to print number of variables
C   --   NOUT   - IN - the output file, <=0 for standard
C   --   NDB    - IN - the database file, <=0 if filename not to be displayed
C   --   TITLE  - IN - the database title
C   --   NUMNP  - IN - the number of nodes
C   --   NDIM   - IN - the number of coordinates per node
C   --   NUMEL  - IN - the number of elements
C   --   NELBLK - IN - the number of element blocks
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   LNPSDF - IN - the length of the node set distribution list
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)
C   --   NVARNS - IN - the number of nodeset variables (if OPTION)
C   --   NVARSS - IN - the number of sideset variables (if OPTION)

      CHARACTER*(*) OPTION
      CHARACTER*80 TITLE

      CHARACTER*132 FILNAM

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         IF (NOUT .GT. 0) WRITE (NOUT, 10000)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         IF (FILNAM(1:1) .NE. ' ') THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10010) FILNAM(:LENSTR(FILNAM))
            ELSE
               WRITE (*, 10010) FILNAM(:LENSTR(FILNAM))
            END IF
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) TITLE(1:LENSTR(TITLE))
         ELSE
            WRITE (*, 10020) TITLE(1:LENSTR(TITLE))
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &         NDIM, NUMNP, NUMEL, NELBLK
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &         NDIM, NUMNP, NUMEL, NELBLK
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         IF (NUMNPS .LE. 0) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10040, IOSTAT=IDUM) NUMNPS
            ELSE
               WRITE (*, 10040, IOSTAT=IDUM) NUMNPS
            END IF
         ELSE
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10040, IOSTAT=IDUM) NUMNPS, LNPSNL, LNPSDF
            ELSE
               WRITE (*, 10040, IOSTAT=IDUM) NUMNPS, LNPSNL, LNPSDF
            END IF
         END IF
         IF (NUMESS .LE. 0) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10050, IOSTAT=IDUM) NUMESS
            ELSE
               WRITE (*, 10050, IOSTAT=IDUM) NUMESS
            END IF
         ELSE
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10050, IOSTAT=IDUM) NUMESS, LESSEL, LESSNL
            ELSE
               WRITE (*, 10050, IOSTAT=IDUM) NUMESS, LESSEL, LESSNL
            END IF
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN
        IF (NOUT .GT. 0) THEN
          WRITE (NOUT, 10060, IOSTAT=IDUM)
     &      NVARGL, NVARNP, NVAREL, NVARNS, NVARSS
        ELSE
          WRITE (*, 10060, IOSTAT=IDUM)
     &      NVARGL, NVARNP, NVAREL, NVARNS, NVARSS
        END IF
      END IF
      
      RETURN

10000  FORMAT (/, 1X, 'DATABASE INITIAL VARIABLES')
10010  FORMAT (/, 1X, 'Database:   ', A)
10020  FORMAT (/, 1X, A)
10030  FORMAT (
     &   /, 1X, 'Number of coordinates per node       =', I10
     &   /, 1X, 'Number of nodes                      =', I10
     &   /, 1X, 'Number of elements                   =', I10
     &   /, 1X, 'Number of element blocks             =', I10
     &   )
10040  FORMAT (
     &   /, 1X, 'Number of node sets                  =', I10, :
     &   /, 1X, '   Length of node list               =', I10,
     &   /, 1X, '   Length of distribution list       =', I10
     &   )
10050  FORMAT
     &   (  1X, 'Number of side sets                  =', I10, :
     &   /, 1X, '   Length of element list            =', I10
     &   /, 1X, '   Length of node list               =', I10
     &   )
10060  FORMAT (
     &   /, 1X, 'Number of global variables           =', I10
     &   /, 1X, 'Number of variables at each node     =', I10
     &   /, 1X, 'Number of variables at each element  =', I10
     &   /, 1X, 'Number of variables at each nodeset  =', I10,
     &   ' (unsupported)'
     &   /, 1X, 'Number of variables at each sideset  =', I10,
     &   ' (unsupported)'
     &   )
      END
