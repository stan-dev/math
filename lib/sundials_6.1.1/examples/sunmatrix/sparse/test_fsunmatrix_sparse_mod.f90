! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2022, Lawrence Livermore National Security
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

module test_fsunmatrix_sparse
  use, intrinsic :: iso_c_binding
  use test_utilities
  implicit none

  integer(C_LONG), parameter :: N = 5

contains

  integer(C_INT) function smoke_tests() result(fails)

    !======== Inclusions ==========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsunmatrix_sparse_mod
    use fnvector_serial_mod

    !======== Declarations ========
    implicit none

    ! local variables
    type(SUNMatrix), pointer :: A, B               ! SUNMatrix
    type(N_Vector),  pointer :: x, y               ! NVectors
    real(C_DOUBLE),  pointer :: matdat(:)          ! matrix data pointer
    integer(C_LONG), pointer :: inddat(:)          ! indices pointer
    integer(C_LONG)          :: lenrw(1), leniw(1) ! matrix real and int work space size

    integer(C_LONG) :: tmp1
    integer(C_INT)  :: tmp2

    fails = 0

    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)

    !===== Calls to interface =====

    ! constructor
    A => FSUNSparseMatrix(N, N, N*N, CSR_MAT, sunctx)
    if (.not. associated(A)) then
      print *,'>>> FAILED - ERROR in FSUNSparseMatrix; halting'
      stop 1
    end if

    ! misc. matrix functions
    tmp1 = FSUNMatGetID_Sparse(A)
    tmp1 = FSUNSparseMatrix_Rows(A)
    tmp1 = FSUNSparseMatrix_Columns(A)
    tmp1 = FSUNSparseMatrix_NNZ(A)
    tmp1 = FSUNSparseMatrix_NP(A)
    tmp2 = FSUNSparseMatrix_SparseType(A)
    matdat => FSUNSparseMatrix_Data(A)
    inddat => FSUNSparseMatrix_IndexValues(A)
    inddat => FSUNSparseMatrix_IndexPointers(A)

    ! matrix operations
    B => FSUNMatClone_Sparse(A)
    if (.not. associated(B)) then
      print *,'>>> FAILED - ERROR in FSUNMatClone_Sparse; halting'
      stop 1
    end if
    fails = fails + FSUNMatZero_Sparse(A)
    fails = fails + FSUNMatCopy_Sparse(A, B)
    fails = fails + FSUNMatScaleAdd_Sparse(ONE, A, B)
    fails = fails + FSUNMatScaleAddI_Sparse(ONE, A)
    fails = fails + FSUNMatMatvec_Sparse(A, x, y)
    fails = fails + FSUNMatSpace_Sparse(A, lenrw, leniw)

    !======= Cleanup ===========
    call FSUNMatDestroy(A)
    call FSUNMatDestroy(B)
    call FN_VDestroy(x)
    call FN_VDestroy(y)

  end function smoke_tests

  integer(C_INT) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsunmatrix_sparse_mod
    use test_sunmatrix

    implicit none

    type(SUNMatrix), pointer :: DA, DI, A, I
    type(N_Vector),  pointer :: x, y
    real(C_DOUBLE),  pointer :: Adata(:), Idata(:), xdata(:), ydata(:)
    integer(C_LONG)          :: ii, jj, tmp1, tmp2

    fails = 0

    ! create dense A and I
    DA => FSUNDenseMatrix(N, N, sunctx)
    DI => FSUNDenseMatrix(N, N, sunctx)

    ! fill A matrix
    Adata => FSUNDenseMatrix_Data(DA)
    do jj=1, N
      do ii=1, N
        Adata((jj-1)*N + ii) = jj*(ii+jj-2)
      end do
    end do

    ! fill identity matrix
    Idata => FSUNDenseMatrix_Data(DI)
    do jj=1, N
      Idata((jj-1)*N + jj) = ONE
    end do

    ! create sparse versions of A and I
    A => FSUNSparseFromDenseMatrix(DA, ZERO, CSR_MAT)
    I => FSUNSparseFromDenseMatrix(DI, ZERO, CSR_MAT)

    ! create vectors
    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)

    ! fill vector x
    xdata => FN_VGetArrayPointer(x)
    do ii=1, N
      xdata(ii) = ONE / ii
    end do

    ! fill vector y
    ydata => FN_VGetArrayPointer(y)
    do ii=1, N
      tmp1 = ii-1
      tmp2 = tmp1 + N - 1
      ydata(ii) = HALF*(tmp2+1-tmp1)*(tmp1+tmp2)
    end do

    fails = fails + Test_FSUNMatGetID(A, SUNMATRIX_SPARSE, 0)
    fails = fails + Test_FSUNMatClone(A, 0)
    fails = fails + Test_FSUNMatCopy(A, 0)
    fails = fails + Test_FSUNMatZero(A, 0)
    fails = fails + Test_FSUNMatScaleAdd(A, I, 0)
    fails = fails + Test_FSUNMatScaleAddI(A, I, 0)
    fails = fails + Test_FSUNMatMatvec(A, x, y, 0)
    fails = fails + Test_FSUNMatSpace(A, 0)

    ! cleanup
    call FSUNMatDestroy(DA)
    call FSUNMatDestroy(DI)
    call FSUNMatDestroy(A)
    call FSUNMatDestroy(I)
    call FN_VDestroy(x)
    call FN_VDestroy(y)

  end function unit_tests

end module

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunmatrix_sparse

  !======== Declarations ========
  implicit none
  integer(C_INT) :: fails = 0

  !============== Introduction =============
  print *, 'Sparse SUNMatrix Fortran 2003 interface test'

  call Test_Init(c_null_ptr)

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed'
    stop 1
  else
    print *, 'SUCCESS: all unit tests passed'
  end if

  call Test_Finalize()

end program main

! exported functions used by test_sunmatrix
integer(C_INT) function check_matrix(A, B, tol) result(fails)
  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fsunmatrix_sparse_mod
  use test_utilities

  implicit none

  type(SUNMatrix)          :: A, B
  real(C_DOUBLE)           :: tol
  real(C_DOUBLE),  pointer :: Adata(:), Bdata(:)
  integer(C_LONG), pointer :: Aidxvals(:), Bidxvals(:)
  integer(C_LONG), pointer :: Aidxptrs(:), Bidxptrs(:)
  integer(C_LONG)          :: i, np, Annz, Bnnz

  fails = 0

  Adata    => FSUNSparseMatrix_Data(A)
  Bdata    => FSUNSparseMatrix_Data(B)
  Aidxvals => FSUNSparseMatrix_IndexValues(A)
  Bidxvals => FSUNSparseMatrix_IndexValues(B)
  Aidxptrs => FSUNSparseMatrix_IndexPointers(A)
  Bidxptrs => FSUNSparseMatrix_IndexPointers(B)

  Annz = FSUNSparseMatrix_NNZ(A)
  Bnnz = FSUNSparseMatrix_NNZ(B)
  np   = FSUNSparseMatrix_NP(A)

  if (FSUNMatGetID(A) /= FSUNMatGetID(B)) then
    fails = 1
    print *, '>>> ERROR: check_matrix: different number of columns'
    return
  end if
  if (FSUNSparseMatrix_SparseType(A) /= FSUNSparseMatrix_SparseType(B)) then
    fails = 1
    print *, '>>> ERROR: check_matrix: different sparse types'
    return
  end if
  if (FSUNSparseMatrix_Rows(A) /= FSUNSparseMatrix_Rows(B)) then
    fails = 1
    print *, '>>> ERROR: check_matrix: different number of rows'
    return
  end if
  if (FSUNSparseMatrix_Columns(A) /= FSUNSparseMatrix_Columns(B)) then
    fails = 1
    print *, '>>> ERROR: check_matrix: different number of columns'
    return
  end if
  if (Annz /= Bnnz) then
    fails = 1
    print *, '>>> ERROR: check_matrix: different number of nonzeros'
    return
  end if

  ! compare sparsity patterns
  do i = 1, np
    if (Aidxptrs(i) /= Bidxptrs(i)) fails = fails + 1
  end do

  if (fails > 0) then
    print *, '>>> ERROR: check_matrix: different indexptrs'
    return
  end if

  do i = 1, Annz
    if (Aidxvals(i) /= Bidxvals(i)) fails = fails + 1
  end do

  if (fails > 0) then
    print *, '>>> ERROR: check_matrix: different indexvals'
    return
  end if

  ! compare matrix values
  do i = 1, Annz
    if (FNEQTOL(Adata(i), Bdata(i), tol) /= 0) fails = fails + 1
  end do

  if (fails > 0) then
    print *, '>>> ERROR: check_matrix: different entries'
    return
  end if


end function check_matrix

integer(C_INT) function check_matrix_entry(A, c, tol) result(fails)
  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fsunmatrix_sparse_mod
  use test_utilities

  implicit none

  type(SUNMatrix)          :: A
  real(C_DOUBLE)           :: c, tol
  real(C_DOUBLE),  pointer :: Adata(:)
  integer(C_LONG), pointer :: Aidxptrs(:)
  integer(C_LONG)          :: i, np

  fails = 0

  Adata    => FSUNSparseMatrix_Data(A)
  Aidxptrs => FSUNSparseMatrix_IndexPointers(A)

  np = FSUNSparseMatrix_NP(A)

  ! compare data
  do i=1, Aidxptrs(np)
    if (FNEQTOL(Adata(i), c, tol) /= 0) then
      fails = fails + 1
      write(*,'(A,I0,A,E14.7,A,E14.7)') &
        'Adata(', i, ') = ', Adata(i), '  c = ', c
    end if
  end do


end function check_matrix_entry

logical function is_square(A) result(res)
  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fsunmatrix_sparse_mod

  implicit none

  type(SUNMatrix) :: A

  if (FSUNSparseMatrix_Rows(A) == FSUNSparseMatrix_Columns(A)) then
    res = .true.
  else
    res = .false.
  end if

end function is_square
