! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
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
    fvec(1) = lamda*yvec(1) + 1.0/(1.0+tn*tn) - lamda*atan(tn);

    ! return success
    ierr = 0
    return

  end function RhsFn

end module ode_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fcvode_mod                    ! Fortran interface to CVODE
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunnonlinsol_fixedpoint_mod  ! Fortran interface to fixed point SUNNonlinearSolver
  use ode_mod                      ! ODE functions

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
  type(c_ptr) :: sunnls        ! sundials fixed-point nonlinear solver
  type(c_ptr) :: cvode_mem     ! CVODE memory

  ! solution vector, neq is set in the ode_functions module
  real(c_double) :: yvec(neq)

  !======= Internals ============

  ! initialize ODE
  tstart = 0.0d0
  tend   = 10.0d0
  tcur   = tstart
  tout   = tstart
  dtout  = 1.0d0
  nout   = ceiling(tend/dtout)

  ! initialize solution vector
  yvec(1) = 0.0d0

  ! create SUNDIALS N_Vector
  sunvec_y = FN_VMake_Serial(neq, yvec)
  if (.not. c_associated(sunvec_y)) print *,'ERROR: sunvec = NULL'

  ! create CVode memory
  cvode_mem = FCVodeCreate(CV_ADAMS)
  if (.not. c_associated(cvode_mem)) print *,'ERROR: cvode_mem = NULL'

  ! initialize CVode
  ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeInit, ierr = ', ierr, '; halting'
     stop
  end if

  ! set relative and absolute tolerances
  rtol = 1.0d-6
  atol = 1.0d-10

  ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
  if (ierr /= 0) then
     write(*,*) 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
     stop
  end if

  ! create fixed point nonlinear solver object
  sunnls = FSUNNonlinSol_FixedPoint(sunvec_y, 0)
  if (.not. c_associated(sunnls)) print *,'ERROR: sunnls = NULL'

  ! attache nonlinear solver object to CVode
  ierr = FCVodeSetNonlinearSolver(cvode_mem, sunnls)
  if (ierr /= 0) then
    write(*,*) 'Error in FCVodeSetNonlinearSolver, ierr = ', ierr, '; halting'
    stop
  end if

  ! Start time stepping
  print *, '   '
  print *, 'Finished initialization, starting time steps'
  print *, '   '
  print *, '       t           y        '
  print *, '----------------------------'
  print '(2x,2(es12.5,1x))', tcur, yvec(1)
  do outstep = 1,nout

     ! call CVode
     tout = min(tout + dtout, tend)
     ierr = FCVode(cvode_mem, tout, sunvec_y, tcur, CV_NORMAL)
     if (ierr /= 0) then
        write(*,*) 'Error in FCVODE, ierr = ', ierr, '; halting'
        stop
     endif

     ! output current solution
     print '(2x,2(es12.5,1x))', tcur, yvec(1)

  enddo

  ! diagnostics output
  call CVodeStats(cvode_mem)
  call FN_VDestroy_Serial(sunvec_y)

  ! clean up
  call FCVodeFree(cvode_mem)

end program main


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
  print *, ' '

  return

end subroutine CVodeStats
