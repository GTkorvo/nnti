!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact Karla Morris  (knmorri@sandia.gov) or
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************

#include "ForTrilinos_config.h"
module FEpetra_SerialComm
  use ForTrilinos_enums ,only : FT_Epetra_Comm_ID,FT_Epetra_SerialComm_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_error ,only : error
  use FEpetra_Comm      ,only : Epetra_Comm
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  !private                     ! Hide everything by default
  !public :: Epetra_SerialComm ! Expose type/constructors/methods

  type ,extends(Epetra_Comm)        :: Epetra_SerialComm 
  contains
    ! !Barrier Methods
    ! procedure         :: barrier
    ! !Broadcast Methods
    ! procedure :: broadcast_double
    ! procedure :: broadcast_int
    ! procedure :: broadcast_long
    ! procedure :: broadcast_char
    ! !Gather Methods
    ! procedure :: gather_double
    ! procedure :: gather_int
    ! procedure :: gather_long
    ! !Sum Methods
    ! procedure :: sum_double
    ! procedure :: sum_int
    ! procedure :: sum_long
    ! !Max/Min Methods
    ! procedure :: max_double
    ! procedure :: max_int
    ! procedure :: max_long
    ! procedure :: min_double
    ! procedure :: min_int
    ! procedure :: min_long
    ! !Parallel Prefix Methods
    ! procedure :: ScanSum_double
    ! procedure :: ScanSum_int
    ! procedure :: ScanSum_long
    ! !Attribute Accessor Methods
    ! procedure :: MyPID
    ! procedure :: NumProc
     !Gather/catter and Directory Constructors
     !I/O methods
  end type

contains
  !> @name Constructor Functions
  !! @{

  !> <BR> Epetra_SerialComm Serial Constructor.  
  type(Epetra_SerialComm) function Epetra_SerialComm()
  end function
 
  !> @name Constructor Functions
  !! @{
  
  !> <BR> Epetra_SerialComm Copy Constructor. 
  type(Epetra_SerialComm) function Epetra_SerialComm(this)
    type(Epetra_SerialComm) ,intent(in) :: this 
  end function

  !> @name Barrier Methods 
  !! @{

  !> <BR> Epetra_SerialComm Barrier function. 
  !> @brief A no-op for a serial communicator.
  !! Implements Epetra_Comm.
  subroutine barrier(this)
    class(Epetra_SerialComm) ,intent(in) :: this
  end subroutine
 
  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_SerialComm Broadcast function. 
  !> @brief A no-op for a serial communicator.
  !!  Implements Epetra_Comm.
  subroutine broadcast(this,MyVals,count,root,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)               ,intent(in)    :: count &
     !< In On entry, contains the length of the list of MyVals. 
    integer(c_int)               ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine
  
  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_SerialComm Broadcast function. 
  !> @brief A no-op for a serial communicator.
  !!  Implements Epetra_Comm.
  subroutine broadcast(this,MyVals,count,root,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)               ,intent(in)    :: count &
     !< In On entry, contains the length of the list of MyVals.
    integer(c_int)               ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information
  end subroutine

  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_SerialComm Broadcast function. 
  !> @brief A no-op for a serial communicator.
  !!  Implements Epetra_Comm.
  subroutine broadcast_long(this,MyVals,count,root,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long),dimension(:) ,intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)               ,intent(in)    :: count &
     !< In On entry, contains the length of the list of MyVals.
    integer(c_int)               ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information
  end subroutine
 
  !> @name Broadcast Methods
  !! @{

  !> <BR> Epetra_SerialComm Broadcast function. 
  !> @brief A no-op for a serial communicator.
  !!  Implements Epetra_Comm.
  subroutine broadcast(this,MyVals,count,root,err)
    class(Epetra_SerialComm)           ,intent(in)    :: this
    character(kind=c_char),dimension(:),intent(inout) :: MyVals &
     !< InOut On entry, the root processor contains the list of values. On exit, all processors will have the same list of values. Note that values must be allocated on all processor before the broadcast.
    integer(c_int)                     ,intent(in)    :: count &
     !< In On entry, contains the length of the list of MyVals.
    integer(c_int)                     ,intent(in)    :: root &
     !< In On entry, contains the processor from which all processors will receive a copy of MyVals.
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
  
  !> @name GatherAll Methods
  !! @{

  !> <BR> Epetra_SerialComm All Gather function. 
  !> @brief A no-op for a serial communicator.
  !!  Implements Epetra_Comm.
 subroutine gather(this,MyVals,AllVals,count,err)
   class(Epetra_SerialComm)     ,intent(in)    :: this
   real(c_double), dimension(:) ,intent(in)    :: MyVals
   real(c_double), dimension(:) ,intent(inout) :: AllVals
   integer(c_int)               ,intent(in)    :: count
   type(error) ,optional, intent(inout) :: err &
   !< Return any error information.
  end subroutine

  subroutine gather_int(this,MyVals,AllVals,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals
    integer(c_int), dimension(:) ,intent(inout) :: AllVals
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  subroutine gather_long(this,MyVals,AllVals,count,err)
    class(Epetra_SerialComm)      ,intent(in)    :: this
    integer(c_long), dimension(:) ,intent(in)    :: MyVals
    integer(c_long), dimension(:) ,intent(inout) :: AllVals
    integer(c_int)                ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  subroutine sum_double(this,PartialSums,GlobalSums,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialSums
    real(c_double), dimension(:) ,intent(inout) :: GlobalSums
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.   
  end subroutine

  subroutine sum_int(this,PartialSums,GlobalSums,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialSums
    integer(c_int), dimension(:) ,intent(inout) :: GlobalSums
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  subroutine sum_long(this,PartialSums,GlobalSums,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialSums
     integer(c_long), dimension(:),intent(inout) :: GlobalSums
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
  
  subroutine max_double(this,PartialMaxs,GlobalMaxs,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMaxs
    real(c_double), dimension(:) ,intent(inout) :: GlobalMaxs
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  subroutine max_int(this,PartialMaxs,GlobalMaxs,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMaxs
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMaxs
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  subroutine max_long(this,PartialMaxs,GlobalMaxs,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMaxs
    integer(c_long), dimension(:),intent(inout) :: GlobalMaxs
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine
  
  subroutine min_double(this,PartialMins,GlobalMins,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: PartialMins
    real(c_double), dimension(:) ,intent(inout) :: GlobalMins
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  subroutine min_int(this,PartialMins,GlobalMins,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: PartialMins
    integer(c_int), dimension(:) ,intent(inout) :: GlobalMins
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  subroutine min_long(this,PartialMins,GlobalMins,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: PartialMins
    integer(c_long), dimension(:),intent(inout) :: GlobalMins
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  subroutine ScanSum_double(this,MyVals,scan_sums,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    real(c_double), dimension(:) ,intent(in)    :: MyVals 
    real(c_double), dimension(:) ,intent(inout) :: scan_sums
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  subroutine ScanSum_int(this,MyVals,scan_sums,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_int), dimension(:) ,intent(in)    :: MyVals 
    integer(c_int), dimension(:) ,intent(inout) :: scan_sums
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
     !< Return any error information.
  end subroutine

  subroutine ScanSum_long(this,MyVals,scan_sums,count,err)
    class(Epetra_SerialComm)     ,intent(in)    :: this
    integer(c_long), dimension(:),intent(in)    :: MyVals 
     integer(c_long), dimension(:),intent(inout) :: scan_sums
    integer(c_int)               ,intent(in)    :: count
    type(error) ,optional, intent(inout) :: err &
    !< Return any error information.
  end subroutine

  integer(c_int) function MyPID(this)
    class(Epetra_SerialComm)     , intent(in) :: this
  end function

  integer(c_int) function NumProc(this)
    class(Epetra_SerialComm)     , intent(in) :: this
  end function

end module 
