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
! fixedpoint SUNNonlinearSolver implementation.
! -----------------------------------------------------------------

module test_fsunnonlinsol_fixedpoint
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use test_utilities

  implicit none

  integer(C_LONG), parameter :: NEQ   = 3      ! number of equations
  integer(C_INT),  parameter :: MAXIT = 10     ! max nonlinear iters.
  real(C_DOUBLE),  parameter :: TOL   = 1.0e-4 ! nonlinear solver tolerance

  real(C_DOUBLE), parameter :: PI = 3.1415926535898

  ! approximate solution
  real(C_DOUBLE) :: Y1 = 0.5d0
  real(C_DOUBLE) :: Y2 = 0.d0
  real(C_DOUBLE) :: Y3 = -PI/6.0d0

contains

  integer(C_INT) function unit_tests() result(retval)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_nonlinearsolver_mod
    use fnvector_serial_mod
    use fsunnonlinsol_fixedpoint_mod
    use fsundials_context_mod

    implicit none

    type(SUNNonlinearSolver), pointer :: NLS         ! test nonlinear solver
    type(N_Vector),           pointer :: x, y0, y, w ! test vectors
    real(C_DOUBLE),           pointer :: ydata(:)
    integer(C_LONG)                   :: niters(1)
    integer(C_INT)                    :: tmp

    x  => FN_VNew_Serial(NEQ, sunctx)
    y0 => FN_VClone(x)
    y  => FN_VClone(x)
    w  => FN_VClone(x)

    ! set weights
    ydata    => FN_VGetArrayPointer(y0)
    ydata(1) = 0.1d0
    ydata(2) = 0.1d0
    ydata(3) = -0.1d0

    call FN_VConst(1.0d0, w)

    ! create and test NLS
    NLS => FSUNNonlinsol_FixedPoint(y, 0, sunctx)

    retval = FSUNNonlinSolSetSysFn(NLS, c_funloc(FPFunction))
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSetSysFn returned ', retval
      return
    end if

    retval = FSUNNonlinSolSetConvTestFn(NLS, c_funloc(ConvTest), c_null_ptr)
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSetConvTestFn returned ', retval
      return
    end if

    retval = FSUNNonlinSolSetMaxIters(NLS, MAXIT)
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSetMaxIters returned ', retval
      return
    end if

    retval = FSUNNonlinSolSolve(NLS, y0, y, w, TOL, 1, c_loc(x))
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSolve returned ', retval
      return
    end if

    ! extract and print solution
    ydata => FN_VGetArrayPointer(y)

    write(*,*) 'Solution:'
    write(*,'(A,E14.7)') 'y1 = ', ydata(1)
    write(*,'(A,E14.7)') 'y2 = ', ydata(2)
    write(*,'(A,E14.7)') 'y3 = ', ydata(3)

    write(*,*) 'Solution Error:'
    write(*,'(A,E14.7)') 'e1 = ', ydata(1) - Y1
    write(*,'(A,E14.7)') 'e2 = ', ydata(2) - Y2
    write(*,'(A,E14.7)') 'e3 = ', ydata(3) - Y3

    retval = FSUNNonlinSolGetNumIters(NLS, niters)
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolGetNumIters returned ', retval
      return
    end if

    write(*,'(A,I0)') 'Number of nonlinear iterations:', niters(1)

    ! cleanup
    call FN_VDestroy(x)
    call FN_VDestroy(y0)
    call FN_VDestroy(y)
    call FN_VDestroy(w)
    tmp = FSUNNonlinSolFree(NLS)

  end function unit_tests

  integer(C_INT) function ConvTest(NLS, y, del, tol, ewt, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_nonlinearsolver_mod

    implicit none

    type(SUNNonlinearSolver) :: NLS
    type(N_Vector)           :: y, del, ewt
    real(C_DOUBLE), value    :: tol
    type(C_PTR), value       :: mem
    real(C_DOUBLE)           :: delnrm

    ! compute the norm of the correction
    delnrm = FN_VMaxNorm(del)

    if (delnrm <= tol) then
      retval = SUN_NLS_SUCCESS  ! converged
    else
      retval = SUN_NLS_CONTINUE ! not converged
    end if

  end function

  ! -----------------------------------------------------------------------------
  ! Nonlinear system
  !
  ! 3x - cos(yz) - 1/2 = 0
  ! x^2 - 81(y+0.1)^2 + sin(z) + 1.06 = 0
  ! exp(-xy) + 20z + (10 pi - 3)/3 = 0
  !
  ! Nonlinear fixed point function
  !
  ! g1(x,y,z) = 1/3 cos(yz) + 1/6
  ! g2(x,y,z) = 1/9 sqrt(x^2 + sin(z) + 1.06) - 0.1
  ! g3(x,y,z) = -1/20 exp(-xy) - (10 pi - 3) / 60
  !
  ! ----------------------------------------------------------------------------
  integer(C_INT) function FPFunction(y, f, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    implicit none

    type(N_Vector)          :: y, f
    type(C_PTR), value      :: mem
    real(C_DOUBLE), pointer :: ydata(:), fdata(:)
    real(C_DOUBLE)          :: y1, y2, y3

    ydata => FN_VGetArrayPointer(y)
    fdata => FN_VGetArrayPointer(f)

    y1 = ydata(1)
    y2 = ydata(2)
    y3 = ydata(3)

    fdata(1) = (1/3.0d0) * cos(y2*y3) + (1/6.0d0)
    fdata(2) = (1/9.0d0) * sqrt(y1*y1 + sin(y3) + 1.06d0) - 0.1d0
    fdata(3) = -(1/20.d0) * exp(-y1*y2) - (10.d0 * PI - 3.0d0) / 60.d0

    retval = 0

  end function

end module

program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunnonlinsol_fixedpoint

  !======== Declarations ========
  implicit none

  integer(C_INT) :: fails = 0

  !============== Introduction =============
  print *, 'fixedpoint SUNNonlinearSolver Fortran 2003 interface test'

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
