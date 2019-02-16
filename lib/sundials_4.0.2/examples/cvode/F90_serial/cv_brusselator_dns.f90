! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! This test simulates a brusselator problem from chemical kinetics.
! This is an ODE system with 3 components, Y = [u,v,w], satisfying
! the equations,
!
!    du/dt = a - (w+1)*u + v*u^2
!    dv/dt = w*u - v*u^2
!    dw/dt = (b-w)/ep - w*u
!
! for t in the interval [0.0, 10.0], with initial conditions
! Y0 = [3.9, 1.1, 2.8], and parameter values a=1.2, b=2.5, and
! ep=1.0e-5
!
! Here, all three components exhibit a rapid transient change during
! the first 0.2 time units, followed by a slow and smooth evolution.
! ------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  ! number of equations
  integer(c_long), parameter :: neq = 3

  ! ODE parameters
  double precision, parameter :: a  = 1.2d0
  double precision, parameter :: b  = 2.5d0
  double precision, parameter :: ep = 1.0d-5

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
       result(ierr) bind(C,name='RhsFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn        ! current time
    type(c_ptr), value    :: sunvec_y  ! solution N_Vector
    type(c_ptr), value    :: sunvec_f  ! rhs N_Vector
    type(c_ptr), value    :: user_data ! user-defined data

    ! pointers to data in SUNDAILS vectors
    real(c_double), pointer :: yvec(:)
    real(c_double), pointer :: fvec(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    call FN_VGetData_Serial(sunvec_y, yvec)
    call FN_VGetData_Serial(sunvec_f, fvec)

    ! fill RHS vector
    fvec(1) = a  -  (yvec(3) + 1.0d0) * yvec(1)  +  yvec(2) * yvec(1) * yvec(1)
    fvec(2) = yvec(3) * yvec(1)  -  yvec(2) * yvec(1) * yvec(1)
    fvec(3) = (b-yvec(3))/ep - yvec(3) * yvec(1)

    ! return success
    ierr = 0
    return

  end function RhsFn

  ! ----------------------------------------------------------------
  ! JacFn: The Jacobian of the ODE hand side function J = df/dy
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, &
       user_data, tmp1, tmp2, tmp3) &
       result(ierr) bind(C,name='JacFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_serial_mod
    use fsunmatrix_dense_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn               ! current time
    type(c_ptr), value    :: sunvec_y         ! current solution N_Vector
    type(c_ptr), value    :: sunvec_f         ! current rhs N_Vector
    type(c_ptr), value    :: sunmat_J         ! Jacobian SUNMatrix
    type(c_ptr), value    :: user_data        ! user-defined data
    type(c_ptr), value    :: tmp1, tmp2, tmp3 ! workspace N_Vectors

    ! pointers to data in SUNDAILS vectors
    real(c_double), pointer :: yvec(:)

    ! pointer to data in SUNDAILS matrix
    real(c_double), pointer :: Jmat(:,:)

    !======= Internals ============

    ! get data array from SUNDIALS vector
    call FN_VGetData_Serial(sunvec_y, yvec)

    ! get data arrays from SUNDIALS vectors
    call FSUNMatGetData_Dense(sunmat_J, Jmat)

    ! fill Jacobian matrix
    Jmat(1,1) = -(yvec(3)+1.0d0) + 2.0d0*yvec(1)*yvec(2);
    Jmat(2,1) = yvec(3) - 2.0d0*yvec(1)*yvec(2);
    Jmat(3,1) = -yvec(3);

    Jmat(1,2) = yvec(1)*yvec(1);
    Jmat(2,2) = -yvec(1)*yvec(1);
    Jmat(3,2) = 0.0d0;

    Jmat(1,3) = -yvec(1);
    Jmat(2,3) = yvec(1);
    Jmat(3,3) = -1.0d0/ep - yvec(1);

    ! return success
    ierr = 0
    return

  end function JacFn

end module ode_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fcvode_mod             ! Fortran interface to CVODE
  use fnvector_serial_mod    ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod      ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod   ! Fortran interface to dense SUNLinearSolver
  use ode_mod               ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: tstart     ! initial time
  real(c_double) :: tend       ! final time
  real(c_double) :: rtol, atol ! relative and absolute tolerance
  real(c_double) :: dtout      ! output time interval
  real(c_double) :: tout       ! output time
  real(c_double) :: tcur       ! current time

  integer(c_int) :: ierr       ! error flag from C functions
  integer(c_int) :: nout       ! number of outputs

  integer :: outstep           ! output loop counter

  type(c_ptr) :: sunvec_y      ! sundials vector
  type(c_ptr) :: sunmat_A      ! sundials matrix
  type(c_ptr) :: sunlinsol_LS  ! sundials linear solver
  type(c_ptr) :: cvode_mem     ! CVODE memory

  ! solution vector, neq is set in the ode_functions module
  real(c_double) :: yvec(neq)

  !======= Internals ============

  ! initialize ODE
  tstart = 0.0d0
  tend   = 10.0d0
  tcur   = tstart
  tout   = tstart
  dtout  = (tend-tstart)/10.d0
  nout   = ceiling(tend/dtout)

  ! initialize solution vector
  yvec(1) = 3.9d0
  yvec(2) = 1.1d0
  yvec(3) = 2.8d0

  ! create SUNDIALS N_Vector
  sunvec_y = FN_VMake_Serial(neq, yvec)
  if (.not. c_associated(sunvec_y)) print *,'ERROR: sunvec = NULL'

  ! create a dense matrix
  sunmat_A = FSUNDenseMatrix(neq, neq);
  if (.not. c_associated(sunmat_A)) print *,'ERROR: sunmat = NULL'

  ! create a dense linear solver
  sunlinsol_LS = FSUNDenseLinearSolver(sunvec_y, sunmat_A);
  if (.not. c_associated(sunlinsol_LS)) print *,'ERROR: sunlinsol = NULL'

  ! create CVode memory
  cvode_mem = FCVodeCreate(CV_BDF)
  if (.not. c_associated(cvode_mem)) print *,'ERROR: cvode_mem = NULL'

  ! initialize CVode
  ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeInit, ierr = ', ierr, '; halting'
     stop
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-5
  atol = 1.0d-10

  ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
     stop
  end if

  ! attach linear solver
  ierr = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A);
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
     stop
  end if

  ! set Jacobian routine
  ierr = FCVodeSetJacFn(cvode_mem, c_funloc(JacFn))
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeSetJacFn, ierr = ', ierr, '; halting'
     stop
  end if

  ! start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '      t           u           v           w'
  print *, '----------------------------------------------------'
  print '(1x,4(es12.5,1x))', tcur, yvec(1), yvec(2), yvec(3)
  do outstep = 1,nout

     ! call CVode
     tout = min(tout + dtout, tend)
     ierr = FCVode(cvode_mem, tout, sunvec_y, tcur, CV_NORMAL)
     if (ierr /= 0) then
        write(*,*) 'Error in FCVODE, ierr = ', ierr, '; halting'
        stop
     endif

     ! output current solution
     print '(1x,4(es12.5,1x))', tcur, yvec(1), yvec(2), yvec(3)

  enddo

  ! diagnostics output
  call CVodeStats(cvode_mem)

  ! clean up
  call FCVodeFree(cvode_mem)
  call FSUNLinSolFree_Dense(sunlinsol_LS)
  call FSUNMatDestroy_Dense(sunmat_A)
  call FN_VDestroy_Serial(sunvec_y)

end program Main


! ----------------------------------------------------------------
! CVodeStats
!
! Print CVODE statstics to stdandard out
! ----------------------------------------------------------------
subroutine CVodeStats(cvode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use fcvode_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: cvode_mem ! solver memory structure

  integer(c_int)  :: ierr ! error flag

  integer(c_long) :: nsteps     ! num steps
  integer(c_long) :: nfevals    ! num function evals
  integer(c_long) :: nlinsetups ! num linear solver setups
  integer(c_long) :: netfails   ! num error test fails

  integer(c_int) :: qlast ! method order in last step
  integer(c_int) :: qcur  ! method order for next step

  real(c_double) :: hinused ! initial step size
  real(c_double) :: hlast   ! last step size
  real(c_double) :: hcur    ! step size for next step
  real(c_double) :: tcur    ! internal time reached

  integer(c_long) :: nniters  ! nonlinear solver iterations
  integer(c_long) :: nncfails ! nonlinear solver fails

  integer(c_long) :: njevals  ! num Jacobian evaluations

  !======= Internals ============

  ! general solver statistics
  ierr = FCVodeGetIntegratorStats(cvode_mem, nsteps, nfevals, nlinsetups, &
       netfails, qlast, qcur, hinused, hlast, hcur, tcur)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeGetIntegratorStats, ierr = ', ierr, '; halting'
     stop
  end if

  ! nonlinear solver statistics
  ierr = FCVodeGetNonlinSolvStats(cvode_mem, nniters, nncfails)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeGetNonlinSolvStats, ierr = ', ierr, '; halting'
     stop
  end if

  ! nonlinear solver statistics
  ierr = FCVodeGetNumJacEvals(cvode_mem, njevals)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeGetNumJacEvals, ierr = ', ierr, '; halting'
     stop
  end if

  print *, ' '
  print *, ' General Solver Stats:'
  print '(4x,A,i9)'    ,'Total internal steps taken =',nsteps
  print '(4x,A,i9)'    ,'Total rhs function calls   =',nfevals
  print '(4x,A,i9)'    ,'Num lin solver setup calls =',nlinsetups
  print '(4x,A,i9)'    ,'Num error test failures    =',netfails
  print '(4x,A,i9)'    ,'Last method order          =',qlast
  print '(4x,A,i9)'    ,'Next method order          =',qcur
  print '(4x,A,es12.5)','First internal step size   =',hinused
  print '(4x,A,es12.5)','Last internal step size    =',hlast
  print '(4x,A,es12.5)','Next internal step size    =',hcur
  print '(4x,A,es12.5)','Current internal time      =',tcur
  print '(4x,A,i9)'    ,'Num nonlinear solver iters =',nniters
  print '(4x,A,i9)'    ,'Num nonlinear solver fails =',nncfails
  print '(4x,A,i9)'    ,'Num Jacobian evaluations   =',njevals
  print *, ' '

  return

end subroutine CVodeStats
