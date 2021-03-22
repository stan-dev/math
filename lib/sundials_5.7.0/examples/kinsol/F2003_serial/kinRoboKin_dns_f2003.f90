! ------------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! This example solves a nonlinear system from robot kinematics.
!
! Source: "Handbook of Test Problems in Local and Global Optimization",
!             C.A. Floudas, P.M. Pardalos et al.
!             Kluwer Academic Publishers, 1999.
! Test problem 6 from Section 14.1, Chapter 14
!
! The nonlinear system is solved by KINSOL using the DENSE linear
! solver.
!
! Constraints are imposed to make all components of the solution
! be within [-1,1].
! ------------------------------------------------------------------

module prob_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  integer(c_long), parameter :: nvar = 8
  integer(c_long), parameter :: neq = 3*nvar
  real(c_double),  parameter :: ftol = 1.d-5
  real(c_double),  parameter :: stol = 1.d-5

contains

  ! ----------------------------------------------------------------
  ! func: The nonlinear system function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function func(sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='func')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_y  ! solution N_Vector
    type(N_Vector)       :: sunvec_f  ! rhs N_Vector
    type(c_ptr),   value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yd(:)
    real(c_double), pointer :: fd(:)

    ! internal variables
    real(c_double) :: x1, x2, x3, x4, x5, x6, x7, x8
    real(c_double) :: l1, l2, l3, l4, l5, l6, l7, l8
    real(c_double) :: u1, u2, u3, u4, u5, u6, u7, u8
    real(c_double) :: eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8
    real(c_double) :: lb1, lb2, lb3, lb4, lb5, lb6, lb7, lb8
    real(c_double) :: ub1, ub2, ub3, ub4, ub5, ub6, ub7, ub8

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yd => FN_VGetArrayPointer(sunvec_y)
    fd => FN_VGetArrayPointer(sunvec_f)

    ! set shortcut variables
    x1 = yd(1)
    x2 = yd(2)
    x3 = yd(3)
    x4 = yd(4)
    x5 = yd(5)
    x6 = yd(6)
    x7 = yd(7)
    x8 = yd(8)
    l1 = yd(9)
    l2 = yd(10)
    l3 = yd(11)
    l4 = yd(12)
    l5 = yd(13)
    l6 = yd(14)
    l7 = yd(15)
    l8 = yd(16)
    u1 = yd(17)
    u2 = yd(18)
    u3 = yd(19)
    u4 = yd(20)
    u5 = yd(21)
    u6 = yd(22)
    u7 = yd(23)
    u8 = yd(24)

    ! Nonlinear equations
    eq1 = - 0.1238d0*x1 + x7 - 0.001637d0*x2 - 0.9338d0*x4 &
          + 0.004731d0*x1*x3 - 0.3578d0*x2*x3 - 0.3571d0
    eq2 = 0.2638d0*x1 - x7 - 0.07745d0*x2 - 0.6734d0*x4 &
        + 0.2238d0*x1*x3 + 0.7623d0*x2*x3 - 0.6022d0
    eq3 = 0.3578d0*x1 + 0.004731d0*x2 + x6*x8
    eq4 = - 0.7623d0*x1 + 0.2238d0*x2 + 0.3461d0
    eq5 = x1*x1 + x2*x2 - 1.d0
    eq6 = x3*x3 + x4*x4 - 1.d0
    eq7 = x5*x5 + x6*x6 - 1.d0
    eq8 = x7*x7 + x8*x8 - 1.d0

    ! Lower bounds ( l_i = 1 + x_i >= 0)
    lb1 = l1 - 1.d0 - x1
    lb2 = l2 - 1.d0 - x2
    lb3 = l3 - 1.d0 - x3
    lb4 = l4 - 1.d0 - x4
    lb5 = l5 - 1.d0 - x5
    lb6 = l6 - 1.d0 - x6
    lb7 = l7 - 1.d0 - x7
    lb8 = l8 - 1.d0 - x8

    ! Upper bounds ( u_i = 1 - x_i >= 0)
    ub1 = u1 - 1.d0 + x1
    ub2 = u2 - 1.d0 + x2
    ub3 = u3 - 1.d0 + x3
    ub4 = u4 - 1.d0 + x4
    ub5 = u5 - 1.d0 + x5
    ub6 = u6 - 1.d0 + x6
    ub7 = u7 - 1.d0 + x7
    ub8 = u8 - 1.d0 + x8

    ! fill in equations
    fd(1) = eq1
    fd(2) = eq2
    fd(3) = eq3
    fd(4) = eq4
    fd(5) = eq5
    fd(6) = eq6
    fd(7) = eq7
    fd(8) = eq8
    fd(9) = lb1
    fd(10) = lb2
    fd(11) = lb3
    fd(12) = lb4
    fd(13) = lb5
    fd(14) = lb6
    fd(15) = lb7
    fd(16) = lb8
    fd(17) = ub1
    fd(18) = ub2
    fd(19) = ub3
    fd(20) = ub4
    fd(21) = ub5
    fd(22) = ub6
    fd(23) = ub7
    fd(24) = ub8

    ! Return success
    ierr = 0
    return

  end function func


  ! ----------------------------------------------------------------
  ! jac: The nonlinear system Jacobian
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function jac(sunvec_y, sunvec_f, sunmat_J, user_data, sunvec_t1, sunvec_t2) &
       result(ierr) bind(C,name='jac')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_y  ! solution N_Vector
    type(N_Vector)       :: sunvec_f  ! rhs N_Vector
    type(SUNMatrix)      :: sunmat_J  ! Jacobian SUNMatrix
    type(c_ptr),   value :: user_data ! user-defined data
    type(N_Vector)       :: sunvec_t1 ! temporary N_Vectors
    type(N_Vector)       :: sunvec_t2

    ! pointers to data in SUNDIALS vector and matrix
    real(c_double), pointer :: yd(:)
    real(c_double), pointer :: J(:,:)

    ! internal variables
    real(c_double) :: x1, x2, x3, x4, x5, x6, x7, x8
    integer :: i

    !======= Internals ============

    ! get data array from SUNDIALS vector and matrix
    yd => FN_VGetArrayPointer(sunvec_y)
    J(1:24, 1:24) => FSUNDenseMatrix_Data(sunmat_J)

    ! set shortcut variables
    x1 = yd(1)
    x2 = yd(2)
    x3 = yd(3)
    x4 = yd(4)
    x5 = yd(5)
    x6 = yd(6)
    x7 = yd(7)
    x8 = yd(8)

    ! --------------------
    ! Nonlinear equations

    ! -0.1238*x1 + x7 - 0.001637*x2 - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571
    J(1,1) = -0.1238d0 + 0.004731d0*x3
    J(1,2) = -0.001637d0 - 0.3578d0*x3
    J(1,3) = 0.004731d0*x1 - 0.3578d0*x2
    J(1,4) = -0.9338d0
    J(1,7) = 1.d0

    ! 0.2638*x1 - x7 - 0.07745*x2 - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022
    J(2,1) = 0.2638d0 + 0.2238d0*x3
    J(2,2) = -0.07745d0 + 0.7623d0*x3
    J(2,3) = 0.2238d0*x1 + 0.7623d0*x2
    J(2,4) = -0.6734d0
    J(2,7) = -1.d0

    ! 0.3578*x1 + 0.004731*x2 + x6*x8
    J(3,1) = 0.3578d0
    J(3,2) = 0.004731d0
    J(3,6) = x8
    J(3,8) = x6

    ! -0.7623*x1 + 0.2238*x2 + 0.3461
    J(4,1) = -0.7623d0
    J(4,2) = 0.2238d0

    ! x1*x1 + x2*x2 - 1
    J(5,1) = 2.d0*x1
    J(5,2) = 2.d0*x2

    ! x3*x3 + x4*x4 - 1
    J(6,3) = 2.d0*x3
    J(6,4) = 2.d0*x4

    ! x5*x5 + x6*x6 - 1
    J(7,5) = 2.d0*x5
    J(7,6) = 2.d0*x6

    ! x7*x7 + x8*x8 - 1
    J(8,7) = 2.d0*x7
    J(8,8) = 2.d0*x8

    ! --------------------
    ! Lower bounds ( l_i = 1 + x_i >= 0)
    do i = 1,8
       J(8+i,i)   = -1.d0
       J(8+i,8+i) = 1.d0
    end do

    ! --------------------
    ! Upper bounds ( u_i = 1 - x_i >= 0)
    do i = 1,8
       J(16+i,i)    = 1.d0
       J(16+i,16+i) = 1.d0
    end do

    ! Return success
    ierr = 0
    return

  end function jac

end module prob_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fkinsol_mod                ! Fortran interface to KINSOL
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use prob_mod                   ! problem-defining functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double)  :: fnormtol, scsteptol
  integer(c_int)  :: ierr
  integer(c_long) :: mset

  type(N_Vector),        pointer :: sunvec_y      ! sundials vectors
  type(N_Vector),        pointer :: sunvec_s
  type(N_Vector),        pointer :: sunvec_c
  type(SUNMatrix),       pointer :: sunmat_J      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver

  type(c_ptr) :: kmem ! KINSOL memory

  ! solution, scaling and constraint vectors; neq is set in the prob_mod module
  real(c_double), dimension(neq) :: y, scale, constraints

  !======= Internals ============

  ! -------------------------
  ! Print problem description

  print *, " "
  print *, "Robot Kinematics Example"
  print *, "8 variables; -1 <= x_i <= 1"
  print *, "KINSOL Problem size: 8 + 2*8 = 24"

  ! -------------------------
  ! Set initial guess, scaling and constraints

  y = 1.d0
  y(1:nvar) = dsqrt(2.d0)/2.d0
  scale = 1.d0
  constraints = 0.d0
  constraints(nvar+1:neq) = 1.d0

  ! -------------------------
  ! Create SUNDIALS vectors for solution, scales, and constraints

  sunvec_y => FN_VMake_Serial(neq, y)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_s => FN_VMake_Serial(neq, scale)
  if (.not. associated(sunvec_s)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_c => FN_VMake_Serial(neq, constraints)
  if (.not. associated(sunvec_c)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  ! -------------------------
  ! Initialize and allocate memory for KINSOL

  kmem = FKINCreate()
  if (.not. c_associated(kmem)) then
     print *, 'ERROR: kmem = NULL'
     stop 1
  end if

  ierr = FKINInit(kmem, c_funloc(func), sunvec_y)
  if (ierr /= 0) then
     print *, 'Error in FKINInit, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set optional inputs

  ierr = FKINSetConstraints(kmem, sunvec_c)
  if (ierr /= 0) then
     print *, 'Error in FKINSetConstraints, ierr = ', ierr, '; halting'
     stop 1
  end if

  fnormtol = ftol
  ierr = FKINSetFuncNormTol(kmem, fnormtol)
  if (ierr /= 0) then
     print *, 'Error in FKINSetFuncNormTol, ierr = ', ierr, '; halting'
     stop 1
  end if

  scsteptol = stol
  ierr = FKINSetScaledStepTol(kmem, scsteptol)
  if (ierr /= 0) then
     print *, 'Error in FKINSetScaledStepTol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Create dense SUNMatrix

  sunmat_J => FSUNDenseMatrix(neq, neq)
  if (.not. associated(sunmat_J)) then
     print *,'ERROR: sunmat = NULL'
     stop 1
  end if

  ! -------------------------
  ! Create dense SUNLinearSolver object

  sunlinsol_LS => FSUNDenseLinearSolver(sunvec_y, sunmat_J)
  if (.not. associated(sunlinsol_LS)) then
     print *,'ERROR: sunlinsol = NULL'
     stop 1
  end if

  ! -------------------------
  ! Attach the matrix and linear solver to KINSOL

  ierr = FKINSetLinearSolver(kmem, sunlinsol_LS, sunmat_J)
  if (ierr /= 0) then
     print *, 'Error in FKINSetLinearSolver, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set the Jacobian function

  ierr = FKINSetJacFn(kmem, c_funloc(jac))
  if (ierr /= 0) then
     print *, 'Error in FKINSetJacFn, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Indicate exact Newton

  mset = 1
  ierr = FKINSetMaxSetupCalls(kmem, mset)
  if (ierr /= 0) then
     print *, 'Error in FKINSetMaxSetupCalls, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Call KINSol to solve problem
  !
  ! arguments: KINSol memory block
  !            Initial guess on input, solution on output
  !            Globalization strategy choice
  !            Scaling vector for the solution
  !            Scaling vector for the residual

  print *, "Initial guess:"
  call PrintOutput(y)

  ierr = FKINSol(kmem, sunvec_y, KIN_LINESEARCH, sunvec_s, sunvec_s)
  if (ierr /= 0) then
     print *, 'Error in FKINSol, ierr = ', ierr, '; halting'
     stop 1
  end if

  print *, " "
  print *, "Computed solution:"
  call PrintOutput(y)


  ! -------------------------
  ! Print final statistics and free memory

  call PrintFinalStats(kmem)
  call FKINFree(kmem)
  ierr = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_J)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_s)
  call FN_VDestroy(sunvec_c)

end program main


! ----------------------------------------------------------------
! PrintOutput: prints solution at selected points
! ----------------------------------------------------------------
subroutine PrintOutput(y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use prob_mod

  !======= Declarations =========
  implicit none

  ! input variable
  real(c_double), dimension(neq) :: y

  ! internal variables
  integer :: i

  !======= Internals ============

  print *, "     l=x+1          x         u=1-x"
  print *, "   ----------------------------------"

  do i = 1,NVAR
     print '(1x,3(f10.6,3x))', y(i+nvar), y(i), y(i+2*nvar)
  end do
  return

end subroutine PrintOutput


! ----------------------------------------------------------------
! PrintFinalStats
!
! Print KINSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(kmem)

  !======= Inclusions ===========
  use iso_c_binding
  use fkinsol_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: kmem

  integer(c_int)  :: ierr
  integer(c_long) :: nni(1), nfe(1), nje(1), nfeD(1)

  !======= Internals ============

  ierr = FKINGetNumNonlinSolvIters(kmem, nni)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumNonlinSolvIters, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumFuncEvals(kmem, nfe)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumFuncEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumJacEvals(kmem, nje)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumJacEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumLinFuncEvals(kmem, nfeD)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumLinFuncEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  print *, ' '
  print *, 'Final Statistics.. '
  print *, ' '
  print '(2(A,i5))'    ,'nni    =', nni,      '    nfe   =', nfe
  print '(2(A,i5))'    ,'nje    =', nje,      '    nfeD  =', nfeD

  return

end subroutine PrintFinalStats
