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
! Newton SUNNonlinearSolver implementation.
! -----------------------------------------------------------------

module test_fsunnonlinsol_newton
  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fsundials_nvector_mod
  use fsundials_linearsolver_mod
  use test_utilities

  implicit none

  integer(C_LONG), parameter :: NEQ   = 3      ! number of equations
  integer(C_INT),  parameter :: MAXIT = 10     ! max nonlinear iters.
  real(C_DOUBLE),  parameter :: TOL   = 1.0e-2 ! nonlinear solver tolerance

  ! approximate solution
  real(C_DOUBLE) :: Y1 = 0.785196933062355226
  real(C_DOUBLE) :: Y2 = 0.496611392944656396
  real(C_DOUBLE) :: Y3 = 0.369922830745872357

  type, private :: IntegratorMem
    type(N_Vector), pointer  :: y0
    type(N_Vector), pointer  :: ycur
    type(N_Vector), pointer  :: ycor
    type(N_Vector), pointer  :: w
    type(N_Vector), pointer  :: x
    type(SUNMatrix), pointer :: A
    type(SUNLinearSolver), pointer :: LS
  end type

contains

  integer(C_INT) function unit_tests() result(retval)
    use, intrinsic :: iso_c_binding
    use fsundials_context_mod
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsundials_linearsolver_mod
    use fsundials_nonlinearsolver_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod
    use fsunlinsol_dense_mod
    use fsunnonlinsol_newton_mod

    implicit none

    type(SUNNonlinearSolver), pointer :: NLS        ! test nonlinear solver
    real(C_DOUBLE),           pointer :: ydata(:)
    integer(C_LONG)                   :: niters(1)
    integer(C_INT)                    :: tmp
    type(IntegratorMem),      pointer :: Imem

    retval = 0

    ! create mock integrator memory
    allocate(Imem)

    ! create vectors
    Imem%y0   => FN_VNew_Serial(NEQ, sunctx)
    Imem%ycur => FN_VClone(Imem%y0)
    Imem%ycor => FN_VClone(Imem%y0)
    Imem%w    => FN_VClone(Imem%y0)
    Imem%x    => FN_VClone(Imem%y0)

    ! set initial guess for the state
    call FN_VConst(HALF, Imem%y0)

    ! set the initial guess for the correction
    call FN_VConst(ZERO, Imem%ycor)

    ! set weights for norm
    call FN_VConst(ONE, Imem%w)

    ! create matrix and linear solver
    Imem%A  => FSUNDenseMatrix(NEQ, NEQ, sunctx)
    Imem%LS => FSUNLinSol_Dense(Imem%y0, Imem%A, sunctx)

    retval = FSUNLinSolInitialize(Imem%LS)
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNLinSolInitialize returned ', retval
      return
    end if

    ! create and test NLS
    NLS => FSUNNonlinsol_Newton(Imem%y0, sunctx)

    retval = FSUNNonlinSolSetSysFn(NLS, c_funloc(Res))
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSetSysFn returned ', retval
      return
    end if

    retval = FSUNNonlinSolSetLSetupFn(NLS, c_funloc(LSetup))
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSetLSetupFn returned ', retval
      return
    end if

    retval = FSUNNonlinSolSetLSolveFn(NLS, c_funloc(LSolve))
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSetLSolveFn returned ', retval
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

    retval = FSUNNonlinSolSolve(NLS, Imem%y0, Imem%ycor, Imem%w, TOL, 1, &
                                c_loc(Imem))
    if (retval /= 0) then
      write(*,'(A,I0)') '   >>> FAIL: FSUNNonlinSolSolve returned ', retval
      return
    end if

    ! update the initial guess with the final correction
    call FN_VLinearSum(ONE, Imem%y0, ONE, Imem%ycor, Imem%ycur)

    ! extract solution data
    ydata => FN_VGetArrayPointer(Imem%ycur)

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
    call FN_VDestroy(Imem%y0)
    call FN_VDestroy(Imem%ycur)
    call FN_VDestroy(Imem%ycor)
    call FN_VDestroy(Imem%w)
    call FN_VDestroy(Imem%x)
    call FSUNMatDestroy(Imem%A)
    tmp = FSUNLinSolFree(Imem%LS)
    tmp = FSUNNonlinSolFree(NLS)
    deallocate(Imem)

  end function unit_tests

  integer(C_INT) function LSetup(jbad, jcur, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod
    use fsundials_nvector_mod

    implicit none

    type(N_Vector), pointer      :: fy, tmp1, tmp2, tmp3
    integer(C_INT), value        :: jbad
    integer(C_INT), dimension(*) :: jcur
    type(C_PTR), value           :: mem
    type(IntegratorMem), pointer :: Imem

    ! set unused parameters to null()
    fy   => null()
    tmp1 => null()
    tmp2 => null()
    tmp3 => null()

    ! get the Integrator memory Fortran type out
    call c_f_pointer(mem, Imem)

    ! compute the Jacobian
    retval = Jac(0.d0, Imem%ycur, fy, Imem%A, C_NULL_PTR, tmp1, tmp2, tmp3)
    if (retval /= 0) return

    ! update Jacobian status
    jcur(1) = 1

    retval = FSUNLinSolSetup(Imem%LS, Imem%A)

  end function

  integer(C_INT) function LSolve(b, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_linearsolver_mod
    use fsundials_nvector_mod

    implicit none

    type(N_Vector)               :: b
    type(C_PTR),         value   :: mem
    type(IntegratorMem), pointer :: Imem

    ! get the Integrator memory Fortran type out
    call c_f_pointer(mem, Imem)

    retval = FSUNLinSolSolve(Imem%LS, Imem%A, Imem%x, b, 0.d0)
    call FN_VScale(1.0d0, Imem%x, b)

  end function

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
    delnrm = FN_VWrmsNorm(del, ewt)

    if (delnrm <= tol) then
      retval = SUN_NLS_SUCCESS  ! converged
    else
      retval = SUN_NLS_CONTINUE ! not converged
    end if

  end function

  integer(C_INT) function Res(ycor, f, mem) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    implicit none

    type(N_Vector)               :: ycor, f
    type(C_PTR), value           :: mem
    real(C_DOUBLE), pointer      :: ydata(:), fdata(:)
    real(C_DOUBLE)               :: y1, y2, y3
    type(IntegratorMem), pointer :: Imem

    ! get the Integrator memory Fortran type out
    call c_f_pointer(mem, Imem)

    ! update the state based on the current correction
    call FN_VLinearSum(ONE, Imem%y0, ONE, ycor, Imem%ycur)

    ydata => FN_VGetArrayPointer(Imem%ycur)
    fdata => FN_VGetArrayPointer(f)

    y1 = ydata(1)
    y2 = ydata(2)
    y3 = ydata(3)

    fdata(1) = y1*y1 + y2*y2 + y3*y3 - 1.0d0
    fdata(2) = 2.0d0 * y1*y1 + y2*y2 - 4.0d0 * y3
    fdata(3) = 3 * y1*y1 - 4.0d0 * y2  + y3*y3

    retval = 0

  end function

  integer(C_INT) function Jac(t, y, fy, J, user_data, tmp1, tmp2, tmp3) &
    result(retval) bind(C)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsunmatrix_dense_mod

    implicit none

    real(C_DOUBLE), value   :: t
    type(N_Vector)          :: y, fy, tmp1, tmp2, tmp3
    type(SUNMatrix)         :: J
    type(C_PTR),    value   :: user_data
    real(C_DOUBLE), pointer :: ydata(:), Jdata(:)
    real(C_DOUBLE)          :: y1, y2, y3

    ydata => FN_VGetArrayPointer(y)
    Jdata => FSUNDenseMatrix_Data(J)

    y1 = ydata(1)
    y2 = ydata(2)
    y3 = ydata(3)

    ! dense matrix has column-major ordering
    Jdata(1:9) = [ TWO*y1, FOUR*y1, SIX*y1, &
                   TWO*y2, TWO*y2, -FOUR, &
                   TWO*y3, -FOUR, TWO*y3 ]

    retval = 0

  end function

end module

program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_fsunnonlinsol_newton

  !======== Declarations ========
  implicit none

  integer(C_INT) :: retval = 0

  !============== Introduction =============
  print *, 'Newton SUNNonlinearSolver Fortran 2003 interface test'

  call Test_Init(c_null_ptr)

  retval = unit_tests()
  if (retval /= 0) then
    print *, 'FAILURE: n unit tests failed'
    stop 1
  else
    print *,'SUCCESS: all unit tests passed'
  end if

  call Test_Finalize()

end program main
