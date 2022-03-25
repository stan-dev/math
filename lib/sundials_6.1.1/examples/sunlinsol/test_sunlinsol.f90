! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! Acknowledgements: These testing routines are based on
! test_sunlinsol.c written by David Gardner @ LLNL and Daniel
! R. Reynolds @ SMU.
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
! These test functions are designed to check the SWIG generated
! Fortran interface to a SUNLinearSolver module implementation.
! -----------------------------------------------------------------

module test_sunlinsol
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fsundials_matrix_mod
  use fsundials_types_mod
  use test_utilities

  implicit none

  ! check_vector routine is provided by implementation specific tests
  integer(C_INT), external :: check_vector

contains

  integer(C_INT) function Test_FSUNLinSolGetType(S, mysunid, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    integer(SUNLinearSolver_Type)  :: mysunid, sunid
    integer(C_INT)                 :: myid

    sunid = FSUNLinSolGetType(S)
    if (sunid /= mysunid) then
      failure = 1
      write(*,*) ">>> FAILED test -- FSUNLinSolGetType, Proc", myid
    else if (myid == 0) then
      failure = 0
      write(*,*) "    PASSED test -- FSUNLinSolGetType"
    end if
  end function Test_FSUNLinSolGetType


  integer(C_INT) function Test_FSUNLinSolLastFlag(S, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    integer(C_INT)                 :: myid
    integer(C_LONG)                :: lastflag

    failure = 0

    ! the only way to fail this test is if the function is NULL,
    ! which will cause a seg-fault
    lastflag = FSUNLinSolLastFlag(S)
    if (myid == 0) then
      write(*,'(A,I0,A)') "     PASSED test -- FSUNLinSolLastFlag (", lastflag, ")"
    end if
  end function Test_FSUNLinSolLastFlag


  integer(C_INT) function Test_FSUNLinSolSpace(S, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    integer(C_INT)                 :: myid
    integer(C_LONG)                :: lenrw(1), leniw(1)

    failure = 0

    ! call FSUNLinSolSpace (failure based on output flag)
    failure = FSUNLinSolSpace(S, lenrw, leniw)
    if (failure /= 0) then
      write(*,*) ">>> FAILED test -- FSUNLinSolSpace, Proc ", myid
    else if (myid == 0) then
      write(*,'(A,I0,A,I0)') "     PASSED test -- FSUNLinSolSpace, lenrw = ", &
        lenrw, " leniw = ",  leniw
    end if

  end function Test_FSUNLinSolSpace


  integer(C_INT) function Test_FSUNLinSolNumIters(S, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    integer(C_INT)                 :: myid
    integer(C_INT)                 :: numiters

    failure = 0

    ! the only way to fail this test is if the function is NULL (segfault will happen)
    numiters = FSUNLinSolNumIters(S)

    if (myid == 0) then
      write(*,'(A,I0,A)') "     PASSED test -- FSUNLinSolNumIters (", numiters, ")"
    end if

  end function Test_FSUNLinSolNumIters


  integer(C_INT) function Test_FSUNLinSolResNorm(S, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    integer(C_INT)                 :: myid
    real(C_DOUBLE)                 :: resnorm

    failure = 0

    resnorm = FSUNLinSolResNorm(S)

    if (resnorm < ZERO) then
      write(*,'(A,E14.7,A,I0)') &
        ">>> FAILED test -- FSUNLinSolSolve returned ", resnorm, ", Proc ", myid
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolResNorm "
    end if

  end function Test_FSUNLinSolResNorm


  integer(C_INT) function Test_FSUNLinSolResid(S, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    integer(C_INT)                 :: myid
    type(N_Vector),        pointer :: resid

    failure = 0

    resid => FSUNLinSolResid(S)

    if (.not. associated(resid)) then
      write(*,*) ">>> FAILED test -- FSUNLinSolResid returned NULL N_Vector, Proc ", myid
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolResid "
    end if

  end function Test_FSUNLinSolResid


  integer(C_INT) function Test_FSUNLinSolSetATimes(S, ATdata, ATimes, myid) &
    result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    type(C_PTR)                    :: ATdata
    type(C_FUNPTR)                 :: ATimes
    integer(C_INT)                 :: myid

    failure = 0

    ! try calling SetATimes routine: should pass/fail based on expected input
    failure = FSUNLinSolSetATimes(S, ATdata, ATimes);

    if (failure /= 0) then
      write(*,'(A,I0,A,I0)') &
        ">>> FAILED test -- FSUNLinSolSetATimes returned ", failure, ", Proc ", myid
      failure = 1
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolSetATimes "
    end if

  end function Test_FSUNLinSolSetATimes


  integer(C_INT) function Test_FSUNLinSolSetPreconditioner(S, Pdata, PSetup, PSolve, myid) &
    result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver), pointer :: S
    type(C_PTR)                    :: Pdata
    type(C_FUNPTR)                 :: PSetup, PSolve
    integer(C_INT)                 :: myid

    ! try calling SetPreconditioner routine: should pass/fail based on expected input
    failure = FSUNLinSolSetPreconditioner(S, Pdata, PSetup, PSolve);

    if (failure /= 0) then
      write(*,'(A,I0,A,I0)') &
        ">>> FAILED test -- FSUNLinSolSetPreconditioner returned ", failure, ", Proc ", myid
      failure = 1
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolSetPreconditioner "
    end if

  end function Test_FSUNLinSolSetPreconditioner


  integer(C_INT) function Test_FSUNLinSolSetScalingVectors(S, s1, s2, myid) &
    result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod
    use fsundials_nvector_mod

    implicit none

    type(SUNLinearSolver) :: S
    type(N_Vector)        :: s1, s2
    integer(C_INT)        :: myid

    failure = 0

    ! try calling SetScalingVectors routine: should pass/fail based on expected input
    failure = FSUNLinSolSetScalingVectors(S, s1, s2)

    if (failure /= 0) then
      write(*,'(A,I0,A,I0)') &
        ">>> FAILED test -- FSUNLinSolSetScalingVectors returned ", failure, ", Proc ", myid
      failure = 1
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolSetScalingVectors "
    end if

  end function Test_FSUNLinSolSetScalingVectors


  integer(C_INT) function Test_FSUNLinSolInitialize(S, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver) :: S
    integer(C_INT)        :: myid

    failure = 0

    failure = FSUNLinSolInitialize(S)

    if (failure /= 0) then
      write(*,'(A,I0,A,I0)') &
        ">>> FAILED test -- FSUNLinSolInitialize returned ", failure, ", Proc ", myid
      failure = 1
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolInitialize "
    end if

  end function Test_FSUNLinSolInitialize

  integer(C_INT) function Test_FSUNLinSolSetup(S, A, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_matrix_mod
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver) :: S
    type(SUNMatrix)       :: A
    integer(C_INT)        :: myid

    failure = 0

    failure = FSUNLinSolSetup(S, A)

    if (failure /= 0) then
      write(*,'(A,I0,A,I0)') &
        ">>> FAILED test -- FSUNLinSolSetup returned ", failure, ", Proc ", myid
      failure = 1
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolSetup "
    end if

  end function Test_FSUNLinSolSetup

  ! ----------------------------------------------------------------------
  ! FSUNLinSolSolve Test
  !
  ! This test must follow Test_FSUNLinSolSetup.  Also, x must be the
  ! solution to the linear system A*x = b (for the original A matrix);
  ! while the 'A' that is supplied to this function should have been
  ! 'setup' by the Test_FSUNLinSolSetup() function prior to this call.
  ! ----------------------------------------------------------------------
  integer(C_INT) function Test_FSUNLinSolSolve(S, A, x, b, tol, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsundials_linearsolver_mod

    implicit none

    type(SUNLinearSolver)   :: S
    type(SUNMatrix)         :: A
    type(N_Vector)          :: x, b
    type(N_Vector), pointer :: y
    real(C_DOUBLE)          :: tol
    integer(C_INT)          :: myid

    failure = 0

    ! clone to create solution vector
    y => FN_VClone(x)
    call FN_VConst(ZERO, y)

    ! perform solve
    failure = FSUNLinSolSolve(S, A, y, b, tol)
    if (failure /= 0) then
      write(*,'(A,I0,A,I0)') &
        ">>> FAILED test -- FSUNLinSolSolve returned ", failure, ", Proc ", myid
      return
    end if

    ! Check solution, and copy y into x for return
    failure = check_vector(x, y, 10.0d0*tol)
    call FN_VScale(ONE, y, x)

    if (failure /= 0) then
      write(*,*) ">>> FAILED test -- FSUNLinSolSolve check, Proc ", myid
    else if (myid == 0) then
      write(*,*) "    PASSED test -- FSUNLinSolSolve"
    end if

    call FN_VDestroy(y)

  end function Test_FSUNLinSolSolve

end module
