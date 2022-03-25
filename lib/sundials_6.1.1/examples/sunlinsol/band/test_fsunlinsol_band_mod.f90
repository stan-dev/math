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
! band SUNLinearSolver implementation.
! -----------------------------------------------------------------

module test_fsunlinsol_band
  use, intrinsic :: iso_c_binding
  use test_utilities
  implicit none

  integer(C_LONG), parameter :: N = 10
  integer(C_LONG), parameter :: mu = 2
  integer(C_LONG), parameter :: ml = 3

contains

  integer(C_INT) function unit_tests() result(fails)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsundials_linearsolver_mod
    use fnvector_serial_mod
    use fsunmatrix_band_mod
    use fsunlinsol_band_mod
    use test_sunlinsol

    implicit none

    type(SUNLinearSolver), pointer :: LS           ! test linear solver
    type(SUNMatrix), pointer :: A                  ! test matrices
    type(N_Vector),  pointer :: x, y, b            ! test vectors
    real(C_DOUBLE),  pointer :: colj(:)            ! matrix column data
    real(C_DOUBLE),  pointer :: xdata(:), Adata(:) ! data arrays
    real(C_DOUBLE)           :: tmpr               ! temporary real value
    integer(C_LONG)          :: j, k
    integer(C_LONG)          :: smu, kstart, kend, offset
    integer(C_INT)           :: tmp

    fails = 0
    smu = 0

    A => FSUNBandMatrix(N, mu, ml, sunctx)
    x => FN_VNew_Serial(N, sunctx)
    y => FN_VNew_Serial(N, sunctx)
    b => FN_VNew_Serial(N, sunctx)

    ! fill A matrix with uniform random data in [0, 1/N)
    Adata => FSUNBandMatrix_Data(A)
    do j=1, N
      offset = smu-mu+1 + j*(smu+ml+1)
      colj(-mu:ml) => Adata(offset:mu+ml)
      kstart = merge(-j, -mu, j < mu)
      kend = merge(N-1-j, ml, j > N-1-ml)
      do k=kstart, kend
        call random_number(tmpr)
        colj(k) = tmpr / N
      end do
    end do

    ! fill x vector with uniform random data in [1, 2)
    xdata => FN_VGetArrayPointer(x)
    do j=1, N
      call random_number(tmpr)
      xdata(j) = ONE + tmpr
    end do

    ! scale/shift matrix to ensure diagonal dominance
    fails = FSUNMatScaleAddI(ONE/(mu+ml+1), A)
    if (fails /= 0) then
      print *, 'FAIL:  FSUNMatScaleAddI failure'
      call FSUNMatDestroy(A)
      call FN_VDestroy(x)
      call FN_VDestroy(y)
      call FN_VDestroy(b)
      return
    end if

    ! create RHS vector for linear solve
    fails = FSUNMatMatvec(A, x, b)
    if (fails /= 0) then
      print *, 'FAIL:  FSUNMatMatvec failure'
      call FSUNMatDestroy(A)
      call FN_VDestroy(x)
      call FN_VDestroy(y)
      call FN_VDestroy(b)
      return
    end if

    ! create band linear solver
    LS => FSUNLinSol_Band(x, A, sunctx)

    ! run tests
    fails = fails + Test_FSUNLinSolInitialize(LS, 0)
    fails = fails + Test_FSUNLinSolSetup(LS, A, 0)
    fails = fails + Test_FSUNLinSolSolve(LS, A, x, b, 100*UNIT_ROUNDOFF, 0)

    fails = fails + Test_FSUNLinSolGetType(LS, SUNLINEARSOLVER_DIRECT, 0)
    fails = fails + Test_FSUNLinSolLastFlag(LS, 0)
    fails = fails + Test_FSUNLinSolSpace(LS, 0)

    ! cleanup
    tmp = FSUNLinSolFree(LS)
    call FSUNMatDestroy(A)
    call FN_VDestroy(x)
    call FN_VDestroy(y)
    call FN_VDestroy(b)

  end function unit_tests

end module

integer(C_INT) function check_vector(X, Y, tol) result(failure)
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use test_utilities

  implicit none
  type(N_Vector)  :: x, y
  real(C_DOUBLE)  :: tol, maxerr
  integer(C_LONG) :: i, xlen, ylen
  real(C_DOUBLE), pointer :: xdata(:), ydata(:)

  failure = 0

  xdata => FN_VGetArrayPointer(x)
  ydata => FN_VGetArrayPointer(y)

  xlen = FN_VGetLength(x)
  ylen = FN_VGetLength(y)

  if (xlen /= ylen) then
    print *, 'FAIL: check_vector: different data array lengths'
    failure = 1
    return
  end if

  do i = 1, xlen
    failure = failure + FNEQTOL(xdata(i), ydata(i), FIVE*tol*abs(xdata(i)))
  end do

  if (failure > 0) then
    maxerr = ZERO
    do i = 1, xlen
      maxerr = max(abs(xdata(i)-ydata(i))/abs(xdata(i)), maxerr)
    end do
    write(*,'(A,E14.7,A,E14.7,A)') &
      "FAIL: check_vector failure: maxerr = ", maxerr, "  (tol = ", FIVE*tol, ")"
  end if

end function check_vector

program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunlinsol_band

  !======== Declarations ========
  implicit none
  integer(C_INT) :: fails = 0

  !============== Introduction =============
  print *, 'Band SUNLinearSolver Fortran 2003 interface test'

  call Test_Init(c_null_ptr)

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed'
    stop 1
  else
    print *,'SUCCESS: all unit tests passed'
  end if

  call Test_Finalize()

end program main
