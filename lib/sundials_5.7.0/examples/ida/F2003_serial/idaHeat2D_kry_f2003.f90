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
! Example problem for IDA: 2D heat equation, serial, GMRES.
!
! This example solves a discretized 2D heat equation problem.
! This version uses the Krylov solver Spgmr.
!
! The DAE system solved is a spatial discretization of the PDE
!          du/dt = d^2u/dx^2 + d^2u/dy^2
! on the unit square. The boundary condition is u = 0 on all edges.
! Initial conditions are given by u = 16 x (1 - x) y (1 - y). The
! PDE is treated with central differences on a uniform M x M grid.
! The values of u at the interior points satisfy ODEs, and
! equations u = 0 at the boundaries are appended, to form a DAE
! system of size N = M^2. Here M = 10.
!
! The system is solved with IDA using the Krylov linear solver
! SPGMR. The preconditioner uses the diagonal elements of the
! Jacobian only. Routines for preconditioning, required by
! SPGMR, are supplied here. The constraints u >= 0 are posed
! for all components. Output is taken at t = 0, .01, .02, .04,
! ..., 10.24. Two cases are run -- with the Gram-Schmidt type
! being Modified in the first case, and Classical in the second.
! The second run uses IDAReInit.
! ------------------------------------------------------------------

module dae_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  integer(c_int),  parameter :: nout  = 11
  integer(c_int),  parameter :: mgrid = 10
  integer(c_long), parameter :: neq   = mgrid*mgrid

  real(c_double) :: dx
  real(c_double) :: coeff
  real(c_double) :: pp(mgrid,mgrid)

contains

  ! ----------------------------------------------------------------
  ! resHeat: The DAE residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function resHeat(tres, sunvec_u, sunvec_up, sunvec_r, user_data) &
       result(ierr) bind(C,name='resHeat')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tres      ! current time
    type(N_Vector)        :: sunvec_u  ! solution N_Vector
    type(N_Vector)        :: sunvec_up ! derivative N_Vector
    type(N_Vector)        :: sunvec_r  ! residual N_Vector
    type(c_ptr),    value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: u(:,:)
    real(c_double), pointer :: up(:,:)
    real(c_double), pointer :: r(:,:)

    ! local variables
    integer(c_long) :: i, j

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    u(1:mgrid, 1:mgrid)  => FN_VGetArrayPointer(sunvec_u)
    up(1:mgrid, 1:mgrid) => FN_VGetArrayPointer(sunvec_up)
    r(1:mgrid, 1:mgrid)  => FN_VGetArrayPointer(sunvec_r)

    ! Initialize r to u, to take care of boundary equations
    r = u

    ! Loop over interior points; set res = up - (central difference)
    do j = 2,mgrid-1
       do i = 2,mgrid-1
          r(i,j) = up(i,j) - coeff*( u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4.d0*u(i,j))
       end do
    end do

    ! return success
    ierr = 0
    return

  end function resHeat

  ! ----------------------------------------------------------------
  ! PsetupHeat: Preconditioner setup routine
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function PSetupHeat(t, sunvec_u, sunvec_up, sunvec_r, cj, prec_data) &
       result(ierr) bind(C,name='PSetupHeat')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    real(c_double), value :: cj        ! step size scaling factor
    type(N_Vector)        :: sunvec_u  ! solution N_Vector
    type(N_Vector)        :: sunvec_up ! derivative N_Vector
    type(N_Vector)        :: sunvec_r  ! residual N_Vector
    type(c_ptr),    value :: prec_data ! preconditioner data

    ! local variables
    real(c_double) :: pelinv

    !======= Internals ============

    ! initialize pp to 1
    pp = 1.d0

    ! Compute the inverse of the preconditioner diagonal elements
    pelinv = 1.d0/(cj + 4.d0*coeff)

    ! set the interior points to the correct value for preconditioning
    pp(2:mgrid-1, 2:mgrid-1) = pelinv

    ! return success
    ierr = 0
    return

  end function PSetupHeat

  ! ----------------------------------------------------------------
  ! PsolveHeat: Preconditioner solve routine
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function PSolveHeat(t, sunvec_u, sunvec_up, sunvec_r, sunvec_rhs, &
       sunvec_sol, cj, delta, prec_data) result(ierr) bind(C,name='PSolveHeat')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t          ! current time
    real(c_double), value :: cj         ! step size scaling factor
    real(c_double), value :: delta      ! desired preconditioner solve accuracy
    type(N_Vector)        :: sunvec_u   ! solution N_Vector
    type(N_Vector)        :: sunvec_up  ! derivative N_Vector
    type(N_Vector)        :: sunvec_r   ! residual N_Vector
    type(N_Vector)        :: sunvec_rhs ! rhs N_Vector
    type(N_Vector)        :: sunvec_sol ! solution N_Vector
    type(c_ptr),    value :: prec_data  ! preconditioner data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: rhs(:,:)
    real(c_double), pointer :: sol(:,:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    rhs(1:mgrid, 1:mgrid) => FN_VGetArrayPointer(sunvec_rhs)
    sol(1:mgrid, 1:mgrid) => FN_VGetArrayPointer(sunvec_sol)

    ! Apply preconditioner to rhs to create sol
    sol = rhs * pp

    ! return success
    ierr = 0
    return

  end function PSolveHeat

end module dae_mod
! ------------------------------------------------------------------


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fida_mod                   ! Fortran interface to IDA
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunlinsol_spgmr_mod       ! Fortran interface to spgmr SUNLinearSolver
  use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use dae_mod                    ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double)  :: rtol, atol, t0, t1, tout, tret(1)
  integer(c_int)  :: retval, iout
  integer(c_long) :: netf(1), ncfn(1), ncfl(1)

  type(N_Vector),        pointer :: sunvec_u     ! sundials solution vector
  type(N_Vector),        pointer :: sunvec_up    ! sundials derivative vector
  type(N_Vector),        pointer :: sunvec_c     ! sundials constraints vector
  type(N_Vector),        pointer :: sunvec_r     ! sundials residual vector
  type(SUNMatrix),       pointer :: sunmat_A     ! sundials matrix (empty)
  type(SUNLinearSolver), pointer :: sunlinsol_LS ! sundials linear solver
  type(c_ptr)                    :: ida_mem      ! IDA memory

  ! solution, residual and constraints vectors, mgrid is set in the dae_mod module
  real(c_double), dimension(mgrid,mgrid) :: uu, up, res, constraints

  !======= Internals ============

  ! Assign parameters in dae_mod
  dx = 1.d0/(mgrid-1)
  coeff = 1.d0/(dx * dx)

  ! create N_Vectors
  sunvec_u => FN_VMake_Serial(neq, uu)
  if (.not. associated(sunvec_u)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_up => FN_VMake_Serial(neq, up)
  if (.not. associated(sunvec_up)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_r => FN_VMake_Serial(neq, res)
  if (.not. associated(sunvec_r)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_c => FN_VMake_Serial(neq, constraints)
  if (.not. associated(sunvec_c)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  ! Initialize solution vectors
  call SetInitialProfile(sunvec_u, sunvec_up, sunvec_r)

  ! Set constraints to all 1's for nonnegative solution values
  constraints = 1.d0

  ! Assign various parameters
  t0   = 0.d0
  t1   = 0.01d0
  rtol = 0.d0
  atol = 1.d-3

  ! Call FIDACreate and FIDAInit to initialize solution
  ida_mem = FIDACreate()
  if (.not. c_associated(ida_mem)) then
     print *, 'ERROR: ida_mem = NULL'
     stop 1
  end if

  retval = FIDASetConstraints(ida_mem, sunvec_c)
  if (retval /= 0) then
     print *, 'Error in FIDASetConstraints, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAInit(ida_mem, c_funloc(resHeat), t0, sunvec_u, sunvec_up)
  if (retval /= 0) then
     print *, 'Error in FIDAInit, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDASStolerances(ida_mem, rtol, atol)
  if (retval /= 0) then
     print *, 'Error in FIDASStolerances, retval = ', retval, '; halting'
     stop 1
  end if

  ! Create the linear solver SUNLinSol_SPGMR with left preconditioning
  ! and the default Krylov dimension
  sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_u, PREC_LEFT, 0)
  if (.not. associated(sunlinsol_LS)) then
     print *, 'ERROR: sunlinsol = NULL'
     stop 1
  end if

  ! IDA recommends allowing up to 5 restarts (default is 0)
  retval = FSUNLinSol_SPGMRSetMaxRestarts(sunlinsol_LS, 5)
  if (retval /= 0) then
     print *, 'Error in FSUNLinSol_SPGMRSetMaxRestarts, retval = ', retval, '; halting'
     stop 1
  end if

  ! Attach the linear solver (will NULL SUNMatrix object)
  sunmat_A => null()
  retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A)
  if (retval /= 0) then
     print *, 'Error in FIDASetLinearSolver, retval = ', retval, '; halting'
     stop 1
  end if

  ! Set the preconditioner solve and setup functions */
  retval = FIDASetPreconditioner(ida_mem, c_funloc(PsetupHeat), c_funloc(PsolveHeat))
  if (retval /= 0) then
     print *, 'Error in FIDASetPreconditioner, retval = ', retval, '; halting'
     stop 1
  end if

  ! Print output heading
  call PrintHeader(rtol, atol)

  !-------------------------------------
  ! CASE I
  !-------------------------------------

  ! Print case number, output table heading, and initial line of table

  print *, " "
  print *, " "
  print *, "Case 1: gsytpe = MODIFIED_GS"
  print *, " "
  print *, "   Output Summary (umax = max-norm of solution)"
  print *, " "
  print *, "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps"
  print *, "----------------------------------------------------------------------"

  ! Loop over output times, call IDASolve, and print results

  tout = t1
  do iout = 1,NOUT
     retval = FIDASolve(ida_mem, tout, tret, sunvec_u, sunvec_up, IDA_NORMAL)
     if (retval < 0) then
        print *, 'Error in FIDASolve, retval = ', retval, '; halting'
        stop 1
     end if
     call PrintOutput(ida_mem, tret(1), uu)
     tout = 2.d0*tout
  end do

  ! Print remaining counters
  retval = FIDAGetNumErrTestFails(ida_mem, netf)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumErrTestFails, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumNonlinSolvConvFails(ida_mem, ncfn)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumLinConvFails(ida_mem, ncfl)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumLinConvFails, retval = ', retval, '; halting'
     stop 1
  end if

  print *, " "
  print '(a,i1)', "Error test failures            = ", netf
  print '(a,i1)', "Nonlinear convergence failures = ", ncfn
  print '(a,i1)', "Linear convergence failures    = ", ncfl

  !-------------------------------------
  ! CASE II
  !-------------------------------------

  ! Re-initialize uu, up

  call SetInitialProfile(sunvec_u, sunvec_up, sunvec_r)

  ! Re-initialize IDA and SPGMR

  retval = FIDAReInit(ida_mem, t0, sunvec_u, sunvec_up)
  if (retval /= 0) then
     print *, 'Error in FIDAReInit, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FSUNLinSol_SPGMRSetGSType(sunlinsol_LS, CLASSICAL_GS)
  if (retval /= 0) then
     print *, 'Error in FSUNLinSol_SPGMRSetGSType, retval = ', retval, '; halting'
     stop 1
  end if

  ! Print case number, output table heading, and initial line of table

  print *, " "
  print *, " "
  print *, "Case 2: gsytpe = CLASSICAL_GS"
  print *, " "
  print *, "   Output Summary (umax = max-norm of solution)"
  print *, " "
  print *, "  time     umax       k  nst  nni  nje   nre   nreLS    h      npe nps"
  print *, "----------------------------------------------------------------------"

  ! Loop over output times, call IDASolve, and print results
  tout = t1
  do iout = 1,NOUT
     retval = FIDASolve(ida_mem, tout, tret, sunvec_u, sunvec_up, IDA_NORMAL)
     if (retval < 0) then
        print *, 'Error in FIDASolve, retval = ', retval, '; halting'
        stop 1
     end if
     call PrintOutput(ida_mem, tret(1), uu)
     tout = 2.d0*tout
  end do

  ! Print remaining counters

  retval = FIDAGetNumErrTestFails(ida_mem, netf)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumErrTestFails, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumNonlinSolvConvFails(ida_mem, ncfn)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumLinConvFails(ida_mem, ncfl)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumLinConvFails, retval = ', retval, '; halting'
     stop 1
  end if

  print *, " "
  print '(a,i1)', "Error test failures            = ", netf
  print '(a,i1)', "Nonlinear convergence failures = ", ncfn
  print '(a,i1)', "Linear convergence failures    = ", ncfl

  ! free memory
  call FIDAFree(ida_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_up)
  call FN_VDestroy(sunvec_r)
  call FN_VDestroy(sunvec_c)

end program main


! ----------------------------------------------------------------
! SetInitialProfile: routine to initialize u and up vectors.
! ----------------------------------------------------------------
subroutine SetInitialProfile(sunvec_u, sunvec_up, sunvec_r)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use dae_mod

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector) :: sunvec_u  ! solution N_Vector
  type(N_Vector) :: sunvec_up ! derivative N_Vector
  type(N_Vector) :: sunvec_r  ! residual N_Vector

  ! pointers to data in SUNDIALS vectors
  real(c_double), pointer :: uu(:,:)
  real(c_double), pointer :: up(:,:)
  real(c_double), pointer :: r(:,:)

  ! local variables
  integer(c_long) :: i, j
  real(c_double)  :: xfact, yfact
  integer(c_int)  :: retval

  !======= Internals ============

  ! get data arrays from SUNDIALS vectors
  uu(1:mgrid, 1:mgrid) => FN_VGetArrayPointer(sunvec_u)
  up(1:mgrid, 1:mgrid) => FN_VGetArrayPointer(sunvec_up)
  r(1:mgrid, 1:mgrid)  => FN_VGetArrayPointer(sunvec_r)

  !======= Internals ============

  ! Initialize uu on all grid points
  do j = 1,mgrid
     yfact = dx * (j-1)
     do i = 1,mgrid
        xfact = dx * (i-1)
        uu(i,j) = 16.d0 * xfact * (1.d0 - xfact) * yfact * (1.d0 - yfact)
     end do
  end do

  ! Initialize up vector to 0
  up = 0.d0

  ! resHeat sets res to negative of ODE RHS values at interior points
  retval = resHeat(0.d0, sunvec_u, sunvec_up, sunvec_r, C_NULL_PTR)

  ! Copy -r into up to get correct interior initial up values
  up = -r

  ! Set up at boundary points to zero
  up(1,:)     = 0.d0
  up(mgrid,:) = 0.d0
  up(:,1)     = 0.d0
  up(:,mgrid) = 0.d0

  return
end subroutine SetInitialProfile


! ----------------------------------------------------------------
! PrintHeader: prints first lines of output (problem description)
! ----------------------------------------------------------------
subroutine PrintHeader(rtol, atol)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use dae_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double) :: rtol, atol

  !======= Internals ============

  print *, " "
  print *, "idaHeat2D_kry: Heat equation, serial example problem for IDA"
  print *, "         Discretized heat equation on 2D unit square."
  print *, "         Zero boundary conditions, polynomial initial conditions."
  print '(2(a,i2),a,i3)', "         Mesh dimensions: ", mgrid, " x ", mgrid, &
       "        Total system size: ", neq
  print *, " "
  print '(2(a,f5.3))', "Tolerance parameters:  rtol = ", rtol,"   atol = ", atol
  print *, "Constraints set to force all solution components >= 0."
  print *, "Linear solver: SPGMR, preconditioner using diagonal elements."

  return
end subroutine PrintHeader


! ----------------------------------------------------------------
! PrintOutput
! ----------------------------------------------------------------
subroutine PrintOutput(ida_mem, t, uu)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fida_mod
  use dae_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  type(c_ptr)    :: ida_mem
  real(c_double) :: t, uu(mgrid,mgrid)

  ! internal variables
  integer(c_int)  :: retval, kused(1)
  integer(c_long) :: nst(1), nni(1), nje(1), nre(1), nreLS(1), nli(1), npe(1), nps(1)
  real(c_double)  :: hused(1), umax

  !======= Internals ============

  umax = maxval(abs(uu))

  retval = FIDAGetLastOrder(ida_mem, kused)
  if (retval /= 0) then
     print *, 'Error in FIDAGetLastOrder, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumSteps(ida_mem, nst)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumNonlinSolvIters(ida_mem, nni)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumNonlinSolvIters, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumResEvals(ida_mem, nre)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumResEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetLastStep(ida_mem, hused)
  if (retval /= 0) then
     print *, 'Error in FIDAGetLastStep, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumJtimesEvals(ida_mem, nje)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumJtimesEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumLinIters(ida_mem, nli)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumLinIters, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumLinResEvals(ida_mem, nreLS)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumLinResEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumPrecEvals(ida_mem, npe)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumPrecEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumPrecSolves(ida_mem, nps)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumPrecSolves, retval = ', retval, '; halting'
     stop 1
  end if


  print '(f5.2,1x,es13.5,1x,i1,2x,3(i3,2x),2(i4,2x),es9.2,2x,2(i3,1x))', &
       t, umax, kused, nst, nni, nje, nre, nreLS, hused(1), npe, nps

end subroutine PrintOutput
