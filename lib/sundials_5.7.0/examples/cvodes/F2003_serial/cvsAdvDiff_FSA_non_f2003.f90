! ------------------------------------------------------------------
! Programmer(s): David J. Gardner, and Cody J. Balos @ LLNL
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
! Based on cvsAdvDiff_FSA_non.c by Scott D. Cohen, Alan C. Hindmarsh,
! George D. Byrne, and Radu Serban @ LLNL
! -----------------------------------------------------------------
! Example problem:
!
! The following is a simple example problem, with the program for
! its solution by CVODES. The problem is the semi-discrete form of
! the advection-diffusion equation in 1-D:
!   du/dt = q1! d^2 u / dx^2 + q2! du/dx
! on the interval 0 <= x <= 2, and the time interval 0 <= t <= 5.
! Homogeneous Dirichlet boundary conditions are posed, and the
! initial condition is:
!   u(x,y,t=0) = x(2-x)exp(2x).
! The PDE is discretized on a uniform grid of size MX+2 with
! central differencing, and with boundary values eliminated,
! leaving an ODE system of size NEQ = MX.
! This program solves the problem with the option for nonstiff
! systems: ADAMS method and functional iteration.
! It uses scalar relative and absolute tolerances.
! Output is printed at t = .5, 1.0, ..., 5.
! Run statistics (optional outputs) are printed at the end.
!
! Optionally, CVODES can compute sensitivities with respect to the
! problem parameters q1 and q2.
! Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
! STAGGERED1) can be used and sensitivities may be included in the
! error test or not (error control set on FULL or PARTIAL,
! respectively).
!
! Execution:
!
! If no sensitivities are desired:
!    % cvsAdvDiff_FSA_non -nosensi
! If sensitivities are to be computed:
!    % cvsAdvDiff_FSA_non -sensi sensi_meth err_con
! where sensi_meth is one of {sim, stg, stg1} and err_con is one of
! {t, f}.
! -----------------------------------------------------------------*/

module ode_problem

  use, intrinsic :: iso_c_binding
  implicit none

  ! problem parameters
  real(c_double),  parameter :: XMAX  = 2.0d0
  real(c_double),  parameter :: T0    = 0.0d0
  real(c_double),  parameter :: T1    = 0.5d0
  real(c_double),  parameter :: DTOUT = 0.5d0
  real(c_double),  parameter :: ATOL  = 1e-5
  integer(c_int),  parameter :: NOUT  = 10
  integer(c_int),  parameter :: NP    = 2
  integer(c_int),  parameter :: NS    = 2
  integer(c_long), parameter :: MX    = 10
  integer(c_long), parameter :: NEQ   = MX

  ! problem constants
  real(c_double) :: ZERO  = 0.d0

  ! problem data
  real(c_double) :: p(0:NP-1)
  real(c_double) :: dx

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
  integer(c_int) function RhsFn(tn, nv_u, nv_udot, user_data) &
    result(ierr) bind(C,name='RhsFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn         ! current time
    type(N_Vector)        :: nv_u       ! solution N_Vector
    type(N_Vector)        :: nv_udot    ! rhs N_Vector
    type(c_ptr),    value :: user_data  ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: u(:)
    real(c_double), pointer :: udot(:)

    ! problem constants
    real(c_double) :: hordc, horac

    ! temp variables
    real(c_double) :: ui, ult, urt, hdiff, hadv
    integer(c_int) :: i

    ! extract needed problem constants
    hordc = p(0)/(dx*dx)
    horac = p(1)/(2.0d0*dx)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    u    => FN_VGetArrayPointer(nv_u)
    udot => FN_VGetArrayPointer(nv_udot)

    ! loop over all grid points
    do i=1, NEQ
      ui = u(i)

      if (i /= 1) then
        ult = u(i-1)
      else
        ult = ZERO
      endif

      if (i /= NEQ) then
        urt = u(i+1)
      else
        urt = ZERO
      endif

      ! set diffusion and avection terms and load into udot
      hdiff   = hordc*(ult - 2.0d0*ui + urt)
      hadv    = horac*(urt - ult)
      udot(i) = hdiff + hadv
    end do

    ! return success
    ierr = 0
    return

  end function RhsFn

  ! Set initial conditions in u vector.
  subroutine SetIC(nv_u)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    implicit none

    type(N_Vector), pointer :: nv_u
    integer(c_int)          :: i
    real(c_double)          :: x
    real(c_double), pointer :: u(:)

    ! Set pointer to data array and get local length of u.
    u => FN_VGetArrayPointer(nv_u)

    ! Load initial profile into u vector
    do i=1, NEQ
      x = i*dx
      u(i) = x*(XMAX - x)*exp(2.0d0*x)
    end do

  end subroutine SetIC

end module ode_problem

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fcvodes_mod                   ! Fortran interface to CVODES
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsundials_nvector_mod         ! Fortran interface to generic N_Vector
  use fsundials_nonlinearsolver_mod ! Fortran interface to generic SUNNonlinearSolver
  use fsunnonlinsol_fixedpoint_mod  ! Fortran interface to fixed point SUNNonlinearSolver
  use ode_problem                   ! ODE defining functions and parameters

  !======= Declarations =========
  implicit none

  ! Local variables
  type(c_ptr)                       :: cvodes_mem
  type(N_Vector),           pointer :: u, uiS
  type(c_ptr)                       :: uS
  type(SUNNonlinearSolver), pointer :: NLS, NLSsens
  integer(c_int)                    :: iout, retval
  real(c_double)                    :: reltol, abstol, tout, t(1)
  integer(c_int)                    :: is
  real(c_double)                    :: pbar(0:NS-1)
  integer(c_int)                    :: plist(0:NS-1)

  ! Command line arguments
  integer(c_int) :: sensi, err_con
  integer(c_int) :: sensi_meth

  ! Process arguments
  call ProcessArgs(sensi, sensi_meth, err_con)

  ! Set problem data
  dx   = XMAX/(MX+1)
  p(0) = 1.0d0
  p(1) = 0.5d0

  ! Allocate and set initial states
  u => FN_VNew_Serial(NEQ)
  if (.not. associated(u)) then
    write(*,*) 'ERROR: FN_VNew_Serial returned NULL'
    stop 1
  endif
  call SetIC(u)

  ! Set integration tolerances
  reltol = ZERO
  abstol = ATOL

  ! Create CVODES object
  cvodes_mem = FCVodeCreate(CV_ADAMS)
  if (.not. c_associated(cvodes_mem)) then
    write(*,*) 'ERROR: cvodes_mem = NULL'
    stop 1
  endif

  ! Initialize CVode
  retval = FCVodeInit(cvodes_mem, c_funloc(RhsFn), T0, u)
  call check_retval(retval, "FCVodeInit")

  ! Set relative and absolute tolerences
  retval = FCVodeSStolerances(cvodes_mem, reltol, abstol)
  call check_retval(retval, "FCVodeSStolerances")

  ! Create fixed point nonlinear solver object
  NLS => FSUNNonlinSol_FixedPoint(u, 0)
  if (.not. associated(NLS)) then
    write(*,*) 'ERROR: FSUNNonlinSol_FixedPoint returned NULL'
    stop 1
  endif

  ! Attach nonlinear solver object to CVode
  retval = FCVodeSetNonlinearSolver(cvodes_mem, NLS)
  call check_retval(retval, "FCVodeSetNonlinearSolver")

  write(*,*) ""
  print '(A,i3)', "1-D advection-diffusion equation, mesh size =", MX

  ! Sensitivity-related settings
  if (sensi /= 0) then

    do is=0, NS-1
      plist(is) = is
      pbar(is)  = p(plist(is))
    end do

    uS = FN_VCloneVectorArray(NS, u)
    if (.not. c_associated(uS)) then
      write(*,*) 'ERROR: FN_VCloneVectorArray returned NULL'
      stop 1
    endif

    do is=0, NS-1
      uiS => FN_VGetVecAtIndexVectorArray(uS, is)
      call FN_VConst(ZERO, uiS)
    end do

    retval = FCVodeSensInit1(cvodes_mem, NS, sensi_meth, c_null_funptr, uS)
    call check_retval(retval, "FCVodeSensInit1")

    retval = FCVodeSensEEtolerances(cvodes_mem)
    call check_retval(retval, "FCVodeSensEEtolerances")

    retval = FCVodeSetSensErrCon(cvodes_mem, err_con)
    call check_retval(retval, "FCVodeSetSensErrCon")

    retval = FCVodeSetSensDQMethod(cvodes_mem, CV_CENTERED, ZERO)
    call check_retval(retval, "FCVodeSetSensDQMethod")

    retval = FCVodeSetSensParams(cvodes_mem, p, pbar, plist)
    call check_retval(retval, "FCVodeSetSensParams")

    ! create sensitivity fixed point nonlinear solver object
    if (sensi_meth == CV_SIMULTANEOUS) then
      NLSsens => FSUNNonlinSol_FixedPointSens(NS+1, u, 0)
    else if (sensi_meth == CV_STAGGERED) then
      NLSsens => FSUNNonlinSol_FixedPointSens(NS, u, 0)
    else
      NLSsens => FSUNNonlinSol_FixedPoint(u, 0)
    endif

    if (.not. associated(NLSsens)) then
      write(*,*) 'ERROR: FSUNNonlinSol_FixedPointSens returned NULL'
      stop 1
    endif

    ! attach nonlinear solver object to CVode
    if (sensi_meth == CV_SIMULTANEOUS) then
      retval = FCVodeSetNonlinearSolverSensSim(cvodes_mem, NLSsens)
    else if (sensi_meth == CV_STAGGERED) then
      retval = FCVodeSetNonlinearSolverSensStg(cvodes_mem, NLSsens)
    else
      retval = FCVodeSetNonlinearSolverSensStg1(cvodes_mem, NLSsens)
    endif

    call check_retval(retval, "FCVodeSetNonlinearSolverSens")

    write(*,'(A)',advance="no") "Sensitivity: YES "
    if (sensi_meth == CV_SIMULTANEOUS) then
      write(*,'(A)',advance="no") "( SIMULTANEOUS +"
    else
      if (sensi_meth == CV_STAGGERED) then
        write(*,'(A)',advance="no") "( STAGGERED +"
      else
        write(*,'(A)',advance="no") "( STAGGERED1 +"
      endif
    endif

    if (err_con /= 0) then
      write(*,'(A)',advance="no") " FULL ERROR CONTROL )"
    else
      write(*,'(A)',advance="no") " PARTIAL ERROR CONTROL )"
    endif

  else

    write(*,'(A)') "Sensitivity: NO "

  endif

  ! In loop over output points, call CVode, print results, test for error

  write(*,*) ""
  write(*,*) ""
  write(*,*) "============================================================"
  write(*,*) "     T     Q       H      NST                    Max norm   "
  write(*,*) "============================================================"

  tout = T1
  do iout = 1, NOUT
    retval = FCVode(cvodes_mem, tout, u, t, CV_NORMAL)
    call check_retval(retval, "FCVode")

    call PrintOutput(cvodes_mem, t, u)

    if (sensi /= 0) then
      retval = FCVodeGetSens(cvodes_mem, t, uS)
      call check_retval(retval, "FCVodeGetSens")
      call PrintOutputS(uS)
    endif

    write(*,*) "------------------------------------------------------------"

    tout = tout + DTOUT
  end do

  ! Print final statistics
  call PrintFinalStats(cvodes_mem, sensi, err_con, sensi_meth)

  ! Free memory
  call FN_VDestroy(u)
  if (sensi /= 0) then
    call FN_VDestroyVectorArray(uS, NS)
  endif

  call FCVodeFree(cvodes_mem)
  retval = FSUNNonlinSolFree(NLS)
  if (sensi /= 0) then
    retval = FSUNNonlinSolFree(NLSsens)
  endif

end program main

! Process and verify arguments.
subroutine ProcessArgs(sensi, sensi_meth, err_con)
  use, intrinsic :: iso_c_binding
  use fcvodes_mod
  implicit none

  integer(c_int)    :: sensi, sensi_meth, err_con
  integer(c_int)    :: argc
  character(len=32) :: arg

  argc       = command_argument_count()
  sensi      = 0
  sensi_meth = -1
  err_con    = 0

  if (argc < 1) then
    call WrongArgs()
  endif

  call get_command_argument(1, arg)
  if (arg == "-nosensi") then
    sensi = 0
  else if (arg == "-sensi") then
    sensi = 1
  else
    call WrongArgs()
  endif

  if (sensi /= 0) then

    if (argc /= 3) then
      call WrongArgs()
    endif

    call get_command_argument(2, arg)
    if (arg == "sim") then
      sensi_meth = CV_SIMULTANEOUS
    else if (arg == "stg") then
      sensi_meth = CV_STAGGERED
    else if (arg == "stg1") then
      sensi_meth = CV_STAGGERED1
    else
      call WrongArgs()
    endif

    call get_command_argument(3, arg)
    if (arg == "t") then
      err_con = 1
    else if (arg == "f") then
      err_con = 0
    else
      call WrongArgs()
    endif

  endif
end subroutine

! Print help.
subroutine WrongArgs()
  write(*,*) ""
  write(*,*) "Usage: ./cvsAdvDiff_FSA_non [-nosensi] [-sensi sensi_meth err_con]"
  write(*,*) "         sensi_meth = sim, stg, or stg1"
  write(*,*) "         err_con    = t or f"
  write(*,*) ""
  call exit(0)
end subroutine

! Print the current number of steps, order, stepsize, and max norm of solution
subroutine PrintOutput(cvodes_mem, t, u)
  use, intrinsic :: iso_c_binding
  use fcvodes_mod
  use fsundials_nvector_mod
  implicit none

  ! subroutine args
  type(c_ptr)    :: cvodes_mem
  real(c_double) :: t(1)
  type(N_Vector) :: u

  ! local variables
  integer(c_int)  :: retval
  integer(c_long) :: nst(1)
  integer(c_int)  :: qu(1)
  real(c_double)  :: hu(1)
  real(c_double)  :: unorm

  retval = FCVodeGetNumSteps(cvodes_mem, nst)
  retval = FCVodeGetLastOrder(cvodes_mem, qu)
  retval = FCVodeGetLastStep(cvodes_mem, hu)

  write(*,'(1x,es9.3,1x,i2,2x,es9.3,i5)') t, qu, hu, nst

  unorm = FN_VMaxNorm(u)
  write(*,'(1x,A,es12.4)') "                                Solution        ", unorm

end subroutine

subroutine PrintOutputS(uS)
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod

  ! subroutine args
  type(c_ptr) :: uS

  ! local variables
  type(N_Vector), pointer :: uiS
  real(c_double)          :: norm

  uiS => FN_VGetVecAtIndexVectorArray(uS, 0)
  norm = FN_VMaxNorm(uiS)
  write(*,'(1x,A,es12.4)') "                                Sensitivity 1   ", norm
  uiS => FN_VGetVecAtIndexVectorArray(uS, 1)
  norm = FN_VMaxNorm(uiS)
  write(*,'(1x,A,es12.4)') "                                Sensitivity 2   ", norm

end subroutine

subroutine PrintFinalStats(cvodes_mem, sensi, err_con, sensi_meth)
  use, intrinsic :: iso_c_binding
  use fcvodes_mod
  implicit none

  ! subroutine args
  type(c_ptr)      :: cvodes_mem
  integer(c_int)   :: sensi, err_con, sensi_meth

  ! local variables
  integer(c_int)  :: retval
  integer(c_long) :: nst(1), nfe(1), nsetups(1), netf(1), nni(1), ncfn(1), &
                     nfSe(1), nfeS(1), nsetupsS(1), netfS(1), nniS(1), ncfnS(1)

  retval = FCVodeGetNumSteps(cvodes_mem, nst)
  retval = FCVodeGetNumRhsEvals(cvodes_mem, nfe)
  retval = FCVodeGetNumLinSolvSetups(cvodes_mem, nsetups)
  retval = FCVodeGetNumErrTestFails(cvodes_mem, netf)
  retval = FCVodeGetNumNonlinSolvIters(cvodes_mem, nni)
  retval = FCVodeGetNumNonlinSolvConvFails(cvodes_mem, ncfn)

  if (sensi /= 0) then

    retval = FCVodeGetSensNumRhsEvals(cvodes_mem, nfSe)
    retval = FCVodeGetNumRhsEvalsSens(cvodes_mem, nfeS)
    retval = FCVodeGetSensNumLinSolvSetups(cvodes_mem, nsetupsS)

    if (err_con /= 0) then
      retval = FCVodeGetSensNumErrTestFails(cvodes_mem, netfS)
    else
      netfS = 0
    endif

    if (sensi_meth == CV_STAGGERED .or. sensi_meth == CV_STAGGERED1) then
      retval = FCVodeGetSensNumNonlinSolvIters(cvodes_mem, nniS)
      retval = FCVodeGetSensNumNonlinSolvConvFails(cvodes_mem, ncfnS)
    else
      nniS  = 0
      ncfnS = 0
    endif

  endif

  write(*,*) ""
  write(*,*) "Final Statistics"
  write(*,*) ""
  write(*,'(1x,A,i9)')        "nst     =", nst
  write(*,'(1x,A,i9)')        "nfe     =", nfe
  write(*,'(1x,A,i9,A,i9)')   "nst     =", netf,  "    nsetups   =", nsetups
  write(*,'(1x,A,i9,A,i9)')   "nni     =", nni,   "    ncfn      =", ncfn

  if (sensi /= 0) then
    write(*,*) ""
    write(*,'(1x,A,i9,A,i9)') "nfSe    =", nfSe,  "    nfeS      =", nfeS
    write(*,'(1x,A,i9,A,i9)') "netfS   =", netfS, "    nsetupsS  =", nsetupsS
    write(*,'(1x,A,i9,A,i9)') "nniS    =", nniS,  "    ncfnS     =", ncfnS
  endif

end subroutine

subroutine check_retval(retval, name)
  use, intrinsic :: iso_c_binding

  character(len=*) :: name
  integer(c_int)   :: retval

  if (retval /= 0) then
    write(*,'(A,A,A)') 'ERROR: ', name,' returned nonzero'
    stop 1
  end if
end subroutine
