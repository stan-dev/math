! ------------------------------------------------------------------
! Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
!                modified by Jean M. Sexton @ LBL
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
! The following is a simple example problem with an analytical
! solution.
!
!   dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
!
! for t in the interval [0.0, 10.0], with initial condition: y=0.
!
! The stiffness of the problem is directly proportional to the
! value of lamda. The value of lamda should be negative to
! result in a well-posed ODE; for values with magnitude larger
! than 100 the problem becomes quite stiff.
! ------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  ! number of equations
  integer(c_long), parameter :: neq = 1

  ! ODE parameters
  double precision, parameter :: lamda = -100.0d0

contains

  ! ----------------------------------------------------------------
  ! RhsFn provides the right hand side function for the
  ! ODE: dy/dt = f(t,y)
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn  ! current time
    type(N_Vector) :: sunvec_y   ! solution N_Vector
    type(N_Vector) :: sunvec_f   ! rhs N_Vector
    type(c_ptr) :: user_data     ! user-defined data

    ! pointers to data in SUNDAILS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yvec => FN_VGetArrayPointer(sunvec_y)
    fvec => FN_VGetArrayPointer(sunvec_f)

    ! fill RHS vector
    fvec(1) = lamda*yvec(1) + 1.0/(1.0+tn*tn) - lamda*atan(tn);

    ! return success
    ierr = 0
    return

  end function RhsFn

end module ode_mod

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                ! Fortran interface to the ARKode module
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep module
  use fsundials_nvector_mod      ! Fortran interface to the generic N_Vector
  use fsundials_matrix_mod       ! Fortran interface to the generic SUNMatrix
  use fsundials_linearsolver_mod ! Fortran interface to the generic SUNLinearSolver
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use ode_mod                    ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: tstart                   ! initial time
  real(c_double) :: tend                     ! final time
  real(c_double) :: rtol, atol               ! relative and absolute tolerance
  real(c_double) :: dtout                    ! output time interval
  real(c_double) :: tout                     ! output time
  real(c_double) :: tcur(1)                  ! current time
  integer(c_int) :: imethod, idefault, pq    ! time step adaptivity parameters
  real(c_double) :: adapt_params(3)          ! time step adaptivity parameters
  integer(c_int) :: ierr                     ! error flag from C functions
  integer(c_int) :: nout                     ! number of outputs
  integer(c_int) :: outstep                  ! output loop counter

  type(N_Vector),        pointer :: sunvec_y    ! sundials vector
  type(SUNMatrix),       pointer :: sunmat_A    ! sundials matrix
  type(SUNLinearSolver), pointer :: sunls       ! sundials linear solver
  type(c_ptr)                    :: arkode_mem  ! ARKODE memory
  real(c_double),        pointer :: yvec(:)     ! underlying vector

  !======= Internals ============

  ! initialize ODE
  tstart = 0.0d0
  tend   = 10.0d0
  tcur   = tstart
  tout   = tstart
  dtout  = 1.0d0
  nout   = ceiling(tend/dtout)

  ! create SUNDIALS N_Vector
  sunvec_y => FN_VNew_Serial(neq)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  yvec => FN_VGetArrayPointer(sunvec_y)

  ! initialize solution vector
  call FN_VConst(0.0d0, sunvec_y)

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_null_funptr, c_funloc(RhsFn), tstart, sunvec_y)
  if (.not. c_associated(arkode_mem)) print *,'ERROR: arkode_mem = NULL'

  ! Tell ARKODE to use a dense linear solver.
  sunmat_A => FSUNDenseMatrix(neq, neq)
  if (.not. associated(sunmat_A)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if

  sunls => FSUNDenseLinearSolver(sunvec_y, sunmat_A)
  if (.not. associated(sunls)) then
     print *, 'ERROR: sunls = NULL'
     stop 1
  end if

  ierr = FARKStepSetLinearSolver(arkode_mem, sunls, sunmat_A)
  if (ierr /= 0) then
    write(*,*) 'Error in FARKStepSetLinearSolver'
    stop 1
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-6
  atol = 1.0d-10

  ierr = FARKStepSStolerances(arkode_mem, rtol, atol)
  if (ierr /= 0) then
     write(*,*) 'Error in FARKStepSStolerances, ierr = ', ierr, '; halting'
     stop 1
  end if

  imethod = 4
  idefault = 1
  pq = 0
  ierr = FARKStepSetAdaptivityMethod(arkode_mem, imethod, idefault, pq, adapt_params)
  if (ierr /= 0) then
     write(*,*) 'Error in FARKStepSetAdaptivityMethod, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '       t           y        '
  print *, '----------------------------'
  print '(2x,2(es12.5,1x))', tcur, yvec(1)
  do outstep = 1,nout

     ! call ARKStep
     tout = min(tout + dtout, tend)
     ierr = FARKStepEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
     if (ierr /= 0) then
        write(*,*) 'Error in FARKODE, ierr = ', ierr, '; halting'
        stop 1
     endif

     ! output current solution
     print '(2x,2(es12.5,1x))', tcur, yvec(1)

  enddo

  ! diagnostics output
  call ARKStepStats(arkode_mem)

  ! clean up
  call FARKStepFree(arkode_mem)
  call FN_VDestroy(sunvec_y)
  call FSUNMatDestroy(sunmat_A)
  ierr = FSUNLinSolFree(sunls)

end program main


! ----------------------------------------------------------------
! ARKStepStats
!
! Print ARKODE statstics to stdandard out
! ----------------------------------------------------------------
subroutine ARKStepStats(arkode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: arkode_mem ! solver memory structure

  integer(c_int)  :: ierr          ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nst_a(1)      ! num steps attempted
  integer(c_long) :: nfe(1)        ! num explicit function evals
  integer(c_long) :: nfi(1)        ! num implicit function evals
  integer(c_long) :: nfevals(1)    ! num function evals
  integer(c_long) :: nlinsetups(1) ! num linear solver setups
  integer(c_long) :: netfails(1)   ! num error test fails

  real(c_double)  :: hinused(1)    ! initial step size
  real(c_double)  :: hlast(1)      ! last step size
  real(c_double)  :: hcur(1)       ! step size for next step
  real(c_double)  :: tcur(1)       ! internal time reached

  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nncfails(1)   ! nonlinear solver fails
  integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations

  !======= Internals ============

  ierr = FARKStepGetNumSteps(arkode_mem, nsteps)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumSteps, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumStepAttempts(arkode_mem, nst_a)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumStepAttempts, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumRhsEvals(arkode_mem, nfe, nfi)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumRhsEvals, retval = ', ierr, '; halting'
     stop 1
  end if
  nfevals=nfe+nfi

  ierr = FARKStepGetActualInitStep(arkode_mem, hinused)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetActualInitStep, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetLastStep(arkode_mem, hlast)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetLastStep, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetCurrentStep(arkode_mem, hcur)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetCurrentStep, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetCurrentTime(arkode_mem, tcur)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetCurrentTime, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumLinSolvSetups(arkode_mem, nlinsetups)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumLinSolvSetups, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumErrTestFails(arkode_mem, netfails)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumErrTestFails, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumNonlinSolvIters(arkode_mem, nniters)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumNonlinSolvIters, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumNonlinSolvConvFails, retval = ', ierr, '; halting'
     stop 1
  end if

  ierr = FARKStepGetNumJacEvals(arkode_mem, njacevals)
  if (ierr /= 0) then
     print *, 'Error in FARKStepGetNumJacEvals, retval = ', ierr, '; halting'
     stop 1
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)'    ,'Total internal steps taken    =',nsteps
  print '(4x,A,i9)'    ,'Total internal steps attempts =',nst_a
  print '(4x,A,i9)'    ,'Total rhs function calls      =',nfevals
  print '(4x,A,i9)'    ,'Num lin solver setup calls    =',nlinsetups
  print '(4x,A,i9)'    ,'Num error test failures       =',netfails
  print '(4x,A,es12.5)','First internal step size      =',hinused
  print '(4x,A,es12.5)','Last internal step size       =',hlast
  print '(4x,A,es12.5)','Next internal step size       =',hcur
  print '(4x,A,es12.5)','Current internal time         =',tcur
  print '(4x,A,i9)'    ,'Num nonlinear solver iters    =',nniters
  print '(4x,A,i9)'    ,'Num nonlinear solver fails    =',nncfails
  print *, ' '

  return

end subroutine ARKStepStats
