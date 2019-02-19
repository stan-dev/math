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
! dense SUNMatrix implementation.
! -----------------------------------------------------------------


program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use fsunmatrix_dense_mod
  use fnvector_serial_mod

  !======== Declarations ========
  implicit none

  ! constants
  real(c_double) :: ONE = 1.d0

  ! local variables
  integer(c_int)  :: fails = 0    ! number of test fails
  integer(c_long) :: N = 100      ! dimensions of SUNMatrix
  type(c_ptr)     :: A, B         ! SUNMatrix
  type(c_ptr)     :: x, y         ! NVectors
  type(c_ptr)     :: matdata      ! Matrix data pointer
  integer(c_long) :: lenrw, leniw ! matrix real and int work space size
  integer(c_long) :: val

  x = FN_VNew_Serial(N)
  y = FN_VNew_Serial(N)

  !======= Introduction =========
  print *,'Dense matrix Fortran 2003 interface test'

  !===== Calls to interface =====
  
  ! constructor
  A = FSUNDenseMatrix(N, N)
  if (.not. c_associated(A)) then
    print *,'>>> FAILED - ERROR in FSUNDenseMatrix; halting'
    stop 1
  end if
  
  ! misc. matrix functions
  val = FSUNMatGetID_Dense(A)
  val = FSUNDenseMatrix_Rows(A)
  val = FSUNDenseMatrix_Columns(A)
  val = FSUNDenseMatrix_LData(A)
  matdata = FSUNDenseMatrix_Data(A)
  matdata = FSUNDenseMatrix_Column(A,N)
  
  ! matrix operations 
  B = FSUNMatClone_Dense(A)
  if (.not. c_associated(B)) then
    print *,'>>> FAILED - ERROR in FSUNMatClone_Dense; halting'
    stop 1
  end if
  fails = fails + FSUNMatZero_Dense(A)
  fails = fails + FSUNMatCopy_Dense(A,B)
  fails = fails + FSUNMatScaleAdd_Dense(ONE, A, B)
  fails = fails + FSUNMatScaleAddI_Dense(ONE, A)
  fails = fails + FSUNMatMatvec_Dense(A, x, y)
  fails = fails + FSUNMatSpace_Dense(A, lenrw, leniw)

  ! destructor
  call FSUNMatDestroy_Dense(A)

  !======= Cleanup ===========
  call FSUNMatDestroy_Dense(B)
  call FN_VDestroy_Serial(x)
  call FN_VDestroy_Serial(y)

  if (fails == 0) then
    print *,'    SUCCESS - all tests passed'
  else
    print *,'    FAILURE - ', fails, ' tests failed'
    stop 1
  end if

end program main

