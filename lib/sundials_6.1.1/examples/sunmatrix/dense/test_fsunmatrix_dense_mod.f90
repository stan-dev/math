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
! dense SUNMatrix implementation.
! -----------------------------------------------------------------

module test_fsunmatrix_dense
  use, intrinsic :: iso_c_binding
  use test_utilities
  implicit none

  integer(C_LONG), parameter :: N = 4

contains

  integer(C_INT) function smoke_tests() result(fails)

    !======== Inclusions ==========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsunmatrix_dense_mod
    use fnvector_serial_mod

    !======== Declarations ========
    implicit none

    ! local variables
    type(SUNMatrix), pointer :: A, B               ! SUNMatrix
    type(N_Vector),  pointer :: x, y               ! NVectors
    real(C_DOUBLE),  pointer :: matdat(:)          ! matrix data pointer
    integer(C_LONG)          :: lenrw(1), leniw(1) ! matrix real and int work space size
    integer(C_LONG)          :: val

    fails = 0

    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)

    !===== Calls to interface =====

    ! constructor
    A => FSUNDenseMatrix(N, N, sunctx)
    if (.not. associated(A)) then
      print *,'>>> FAILED - ERROR in FSUNDenseMatrix; halting'
      stop 1
    end if

    ! misc. matrix functions
    val = FSUNMatGetID_Dense(A)
    val = FSUNDenseMatrix_Rows(A)
    val = FSUNDenseMatrix_Columns(A)
    val = FSUNDenseMatrix_LData(A)
    matdat => FSUNDenseMatrix_Data(A)
    matdat => FSUNDenseMatrix_Column(A,N)

    ! matrix operations
    B => FSUNMatClone_Dense(A)
    if (.not. associated(B)) then
      print *,'>>> FAILED - ERROR in FSUNMatClone_Dense; halting'
      stop 1
    end if
    fails = fails + FSUNMatZero_Dense(A)
    fails = fails + FSUNMatCopy_Dense(A,B)
    fails = fails + FSUNMatScaleAdd_Dense(ONE, A, B)
    fails = fails + FSUNMatScaleAddI_Dense(ONE, A)
    fails = fails + FSUNMatMatvec_Dense(A, x, y)
    fails = fails + FSUNMatSpace_Dense(A, lenrw, leniw)

    !======= Cleanup ===========
    call FSUNMatDestroy_Dense(A)
    call FSUNMatDestroy_Dense(B)
    call FN_VDestroy_Serial(x)
    call FN_VDestroy_Serial(y)

  end function smoke_tests

  integer(C_INT) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use test_sunmatrix

    implicit none

    type(SUNMatrix), pointer :: A, I
    type(N_Vector),  pointer :: x, y
    real(C_DOUBLE),  pointer :: Adata(:), Idata(:), xdata(:), ydata(:)
    integer(C_LONG)          :: ii, jj, tmp1, tmp2

    fails = 0

    A => FSUNDenseMatrix(N, N, sunctx)
    I => FSUNDenseMatrix(N, N, sunctx)
    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)

    ! fill matrix A
    Adata => FSUNDenseMatrix_Data(A)
    do jj=1, N
      do ii=1, N
        Adata((jj-1)*N + ii) = jj*(ii+jj-2)
      end do
    end do

    ! fill matrix I (identity)
    Idata => FSUNDenseMatrix_Data(I)
    do jj=1, N
        Idata((jj-1)*N + jj) = ONE
    end do

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

    fails = fails + Test_FSUNMatGetID(A, SUNMATRIX_DENSE, 0)
    fails = fails + Test_FSUNMatClone(A, 0)
    fails = fails + Test_FSUNMatCopy(A, 0)
    fails = fails + Test_FSUNMatZero(A, 0)
    fails = fails + Test_FSUNMatScaleAdd(A, I, 0)
    fails = fails + Test_FSUNMatScaleAddI(A, I, 0)
    fails = fails + Test_FSUNMatMatvec(A, x, y, 0)
    fails = fails + Test_FSUNMatSpace(A, 0)

    ! cleanup
    call FSUNMatDestroy(A)
    call FSUNMatDestroy(I)
    call FN_VDestroy(x)
    call FN_VDestroy(y)

  end function unit_tests

end module

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunmatrix_dense

  !======== Declarations ========
  implicit none
  integer(C_INT) :: fails = 0

  !============== Introduction =============
  print *, 'Dense SUNMatrix Fortran 2003 interface test'

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
  use fsunmatrix_dense_mod
  use test_utilities

  implicit none

  type(SUNMatrix) :: A, B
  real(C_DOUBLE)  :: tol
  real(C_DOUBLE), pointer :: Adata(:), Bdata(:)
  integer(C_LONG) :: Aldata, Bldata, i

  fails = 0

  ! get data pointers
  Adata => FSUNDenseMatrix_Data(A)
  Bdata => FSUNDenseMatrix_Data(B)

  ! get and check data lengths
  Aldata = FSUNDenseMatrix_LData(A)
  Bldata = FSUNDenseMatrix_LData(B)

  if (Aldata /= Bldata) then
    print *, ">>> ERROR: check_matrix: Different data array lengths"
    fails = 1
    return
  end if

  ! compare data
  do i=1, Aldata
    fails = fails + FNEQTOL(Adata(i), Bdata(i), tol)
  end do

end function check_matrix

integer(C_INT) function check_matrix_entry(A, c, tol) result(fails)
  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fsunmatrix_dense_mod
  use test_utilities

  implicit none

  type(SUNMatrix) :: A
  real(C_DOUBLE)  :: c, tol
  real(C_DOUBLE), pointer :: Adata(:)
  integer(C_LONG) :: Aldata, i

  fails = 0

  ! get data pointers
  Adata => FSUNDenseMatrix_Data(A)

  ! get and check data lengths
  Aldata = FSUNDenseMatrix_LData(A)

  ! compare data
  do i=1, Aldata
    fails = fails + FNEQTOL(Adata(i), c, tol)
  end do

  if (fails > ZERO) then
    print *, ">>> ERROR: check_matrix_entry failures: "
    do i=1, Aldata
      if (FNEQTOL(Adata(i), c, tol) /= 0) then
        write(*,'(A,I0,A,E14.7,A,E14.7)') &
          "Adata[ ", i, "] =", Adata(i) ," c = ", c
      end if
    end do
  end if

end function check_matrix_entry

logical function is_square(A) result(res)
  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fsunmatrix_dense_mod

  implicit none

  type(SUNMatrix) :: A

  if (FSUNDenseMatrix_Rows(A) == FSUNDenseMatrix_Columns(A)) then
    res = .true.
  else
    res = .false.
  end if

end function is_square
