! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------
! This file tests the Fortran 2003 interface to the SUNDIALS 
! banded SUNMatrix implementation.
! -----------------------------------------------------------------


program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use fsunmatrix_band_mod
  use fnvector_serial_mod

  !======== Declarations ========
  implicit none

  ! constants
  real(c_double) :: ONE  = 1.d0

  ! local variables
  integer(c_int)  :: fails = 0    ! number of test fails
  integer(c_long) :: N            ! dimensions of SUNMatrix
  integer(c_long) :: uband, lband ! bandwidth of band SUNMatrix
  type(c_ptr)     :: A, B         ! Band SUNMatrix
  type(c_ptr)     :: x, y         ! NVectors
  type(c_ptr)     :: matdata      ! Matrix data pointer
  integer(c_long) :: lenrw, leniw ! matrix real and int work space size
  integer(c_long) :: val

  N = 100
  uband = 0
  lband = 0

  x = FN_VNew_Serial(N)
  y = FN_VNew_Serial(N)
  
  !======= Introduction =========
  print *,'Band matrix Fortran 2003 interface test'

  !===== Calls to interface =====
  
  ! constructors
  A = FSUNBandMatrix(N, uband, lband, uband)
  if (.not. c_associated(A)) then
    print *,'>>> FAILED - ERROR in FSUNBandMatrix; halting'
    stop 1
  end if
 
  ! misc. matrix functions
  val = FSUNBandMatrix_Rows(A)
  val = FSUNBandMatrix_Columns(A)
  val = FSUNBandMatrix_LowerBandwidth(A)
  val = FSUNBandMatrix_UpperBandwidth(A)
  val = FSUNBandMatrix_StoredUpperBandwidth(A)
  val = FSUNBandMatrix_LDim(A)
  matdata = FSUNBandMatrix_Data(A)
  matdata = FSUNBandMatrix_Column(A, N)

  ! matrix operations 
  B = FSUNMatClone_Band(A)
  if (.not. c_associated(B)) then
    print *,'>>> FAILED - ERROR in FSUNMatClone_Band; halting'
    stop 1
  end if
  val = FSUNMatGetID_Band(A)
  fails = fails + FSUNMatZero_Band(A)
  fails = fails + FSUNMatCopy_Band(A,B)
  fails = fails + FSUNMatScaleAdd_Band(ONE, A, B)
  fails = fails + FSUNMatScaleAddI_Band(ONE, A)
  fails = fails + FSUNMatMatvec_Band(A, x, y)
  fails = fails + FSUNMatSpace_Band(A, lenrw, leniw)

  ! destructor
  call FSUNMatDestroy_Band(A)

  !======= Cleanup ===========
  call FSUNMatDestroy_Band(B)
  call FN_VDestroy_Serial(x)
  call FN_VDestroy_Serial(y)

  if (fails == 0) then
    print *,'    SUCCESS - all tests passed'
  else
    print *,'    FAILURE - ', fails, ' tests failed'
    stop 1
  end if

end program main

