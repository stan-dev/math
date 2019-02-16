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
! sparse SUNMatrix implementation.
! -----------------------------------------------------------------


program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use fsunmatrix_sparse_mod
  use fsunmatrix_dense_mod
  use fsunmatrix_band_mod
  use fnvector_serial_mod

  !======== Declarations ========
  implicit none

  ! constants
  real(c_double) :: ZERO = 0.d0
  real(c_double) :: ONE  = 1.d0

  ! local variables
  integer(c_int)  :: fails = 0    ! number of test fails
  integer(c_int)  :: typ          ! type of Sparse matrix
  integer(c_long) :: N            ! dimensions of SUNMatrix
  integer(c_long) :: uband, lband ! bandwidth of band SUNMatrix
  type(c_ptr)     :: A, B         ! Sparse SUNMatrix
  type(c_ptr)     :: D, E         ! Dense and Banded SUNMatrix
  type(c_ptr)     :: x, y         ! NVectors
  type(c_ptr)     :: matdata      ! Matrix data pointer
  integer(c_long) :: lenrw, leniw ! matrix real and int work space size
  integer(c_long) :: val

  N = 100
  uband = 0
  lband = 0

  x = FN_VNew_Serial(N)
  y = FN_VNew_Serial(N)
  
  D = FSUNDenseMatrix(N, N)
  E = FSUNBandMatrix(N, uband, lband, uband)

  !======= Introduction =========
  print *,'Sparse matrix Fortran 2003 interface test'

  !===== Calls to interface =====
  
  ! constructors
  A = FSUNSparseMatrix(N, N, N, CSC_MAT)
  if (.not. c_associated(A)) then
    print *,'>>> FAILED - ERROR in FSUNSparseMatrix; halting'
    stop 1
  end if
  call FSUNMatDestroy_Sparse(A)
  
  A = FSUNSparseFromDenseMatrix(D, ZERO, CSR_MAT)
  if (.not. c_associated(A)) then
    print *,'>>> FAILED - ERROR in FSUNSparseFromDenseMatrix; halting'
    stop 1
  end if
  call FSUNMatDestroy_Sparse(A)
  
  A = FSUNSparseFromBandMatrix(E, ZERO, CSC_MAT)
  if (.not. c_associated(A)) then
    print *,'>>> FAILED - ERROR in FSUNSparseFromBandMatrix; halting'
    stop 1
  end if
  
  ! misc. matrix functions
  val = FSUNSparseMatrix_Rows(A)
  val = FSUNSparseMatrix_Columns(A)
  val = FSUNSparseMatrix_NNZ(A)
  val = FSUNSparseMatrix_NP(A)
  typ = FSUNSparseMatrix_SparseType(A)
  matdata = FSUNSparseMatrix_Data(A)
  matdata = FSUNSparseMatrix_IndexValues(A)
  matdata = FSUNSparseMatrix_IndexPointers(A)
  fails = fails + FSUNSparseMatrix_Realloc(A)
  fails = fails + FSUNSparseMatrix_Reallocate(A, N)
  
  ! matrix operations 
  B = FSUNMatClone_Sparse(A)
  if (.not. c_associated(B)) then
    print *,'>>> FAILED - ERROR in FSUNMatClone_Sparse; halting'
    stop 1
  end if
  val = FSUNMatGetID_Sparse(A)
  fails = fails + FSUNMatZero_Sparse(A)
  fails = fails + FSUNMatCopy_Sparse(A,B)
  fails = fails + FSUNMatScaleAdd_Sparse(ONE, A, B)
  fails = fails + FSUNMatScaleAddI_Sparse(ONE, A)
  fails = fails + FSUNMatMatvec_Sparse(A, x, y)
  fails = fails + FSUNMatSpace_Sparse(A, lenrw, leniw)

  ! destructor
  call FSUNMatDestroy_Sparse(A)

  !======= Cleanup ===========
  call FSUNMatDestroy_Sparse(B)
  call FSUNMatDestroy_Dense(D)
  call FSUNMatDestroy_Band(E)
  call FN_VDestroy_Serial(x)
  call FN_VDestroy_Serial(y)

  if (fails == 0) then
    print *,'    SUCCESS - all tests passed'
  else
    print *,'    FAILURE - ', fails, ' tests failed'
    stop 1
  end if

end program main

