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

C=======================================================================
      SUBROUTINE RELBLK (IELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, *)
C=======================================================================

C   $Id: relblk.f,v 1.1 1990/08/20 12:22:44 gdsjaar Exp $
C   $Log: relblk.f,v $
C   Revision 1.1  1990/08/20 12:22:44  gdsjaar
C   Initial revision
C

C   --*** RELBLK *** (GEN3D) Read database element block misc.
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RELBLK reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   IELB - IN - the element block number (for errors)
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - OUT - the element connectivity for this block (4 nodes always)
C   --   ATRIB - OUT - the attributes for this block (packed)
C   --   * - return statement if end of file or read error
C   --
C   --Common Variables:
C   --   Uses NDBIN of /DBASE/
C   --
C   --Database must be positioned at start of element block misc. information
C   --upon entry; upon exit at end of element block misc. information.

      INCLUDE 'gs_dbase.blk'

      INTEGER LINK(4,NUMELB)
      REAL ATRIB(*)

      CHARACTER*5 STRA

      READ (NDBIN, IOSTAT=IERR, END=10, ERR=10)
     &   ((LINK(I,NE), I=1,NUMLNK), NE=1,NUMELB)
      READ (NDBIN, IOSTAT=IERR, END=10, ERR=10)
     &   (ATRIB(I), I=1,NUMATR*NUMELB)

      RETURN

   10 CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL', 'Reading BLOCK ' // STRA(:LSTRA))
      RETURN 1
      END
