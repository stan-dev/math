! -------------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! Based on ark_kpr_mri.c in examples/arkode/C_serial.
! -------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2022, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -------------------------------------------------------------------
! Multirate nonlinear Kvaerno-Prothero-Robinson ODE test problem:
!
!    [u]' = [ G  e ] [(-1+u^2-r)/(2u)] + [      r'(t)/(2u)        ]
!    [v]    [ e -1 ] [(-2+v^2-s)/(2v)]   [ s'(t)/(2*sqrt(2+s(t))) ]
!         = [ fs(t,u,v) ]
!           [ ff(t,u,v) ]
!
! where r(t) = 0.5*cos(t),  s(t) = cos(w*t),  0 < t < 5.
!
! This problem has analytical solution given by
!    u(t) = sqrt(1+r(t)),  v(t) = sqrt(2+s(t)).
!
! We use the parameters:
!   e = 0.5 (fast/slow coupling strength) [default]
!   G = -1e2 (stiffness at slow time scale) [default]
!   w = 100  (time-scale separation factor) [default]
!   hs = 0.01 (slow step size) [default]
!
! The stiffness of the slow time scale is essentially determined
! by G, for |G| > 50 it is 'stiff' and ideally suited to a
! multirate method that is implicit at the slow time scale.
!
! We select the MRI method to use based on an additional input,
! solve_type with options (slow type-order/fast type-order):
! 0. exp-3/exp-3 (standard MIS) [default]
! 1. none/exp-3 (no slow, explicit fast)
! 2. none/dirk-3 (no slow, dirk fast)
! 3. exp-3/none (explicit slow, no fast)
! 4. dirk-2/none (dirk slow, no fast) -- solve-decoupled
! 5. exp-4/exp-4 (MRI-GARK-ERK45a / ERK-4-4)
! 6. exp-4/exp-3 (MRI-GARK-ERK45a / ERK-3-3)
! 7. dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled
! 8. ars343/exp-3 (IMEX-MRI3b / ERK-3-3) -- solve decoupled
! 9. imexark4/exp-4 (IMEX-MRI4/ ERK-4-4) -- solve decoupled
!
! We note that once we have methods that are IMEX at the slow time
! scale, the nonstiff slow term,  [ r'(t)/(2u) ], can be treated
! explicitly.
!
! The program should be run with arguments in the following order:
!   $ a.out solve_type h G w e
! Not all arguments are required, but these must be omitted from
! end-to-beginning, i.e. any one of
!   $ a.out solve_type h G w
!   $ a.out solve_type h G
!   $ a.out solve_type h
!   $ a.out solve_type
!   $ a.out
! are acceptable.  We require:
!   * 0 <= solve_type <= 9
!   * 0 < h < 1/|G|
!   * G < 0.0
!   * w >= 1.0
!
! This program solves the problem with the MRI stepper. Outputs are
! printed at equal intervals of 0.1 and run statistics are printed
! at the end.
!--------------------------------------------------------------------

module ode_mod

  use, intrinsic :: iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod
  use farkode_mristep_mod
  use fnvector_serial_mod
  use fsunmatrix_dense_mod
  use fsunlinsol_dense_mod
  use fsundials_context_mod
  use fsundials_matrix_mod
  use fsundials_linearsolver_mod
  use fsundials_nvector_mod

  implicit none

  ! Constants
  real(c_double) :: ZERO = 0.0d0
  real(c_double) :: ONE  = 1.0d0
  real(c_double) :: TWO  = 2.0d0

  ! general problem parameters
  real(c_double), parameter  :: T0          = 0.0d0    ! initial time
  real(c_double), parameter  :: Tf          = 5.0d0    ! final time
  real(c_double), parameter  :: dTout       = 0.1d0    ! time between outputs
  integer(c_long), parameter :: NEQ         = 2        ! number of dependent vars.
  integer(c_int), parameter  :: Nt          = ceiling(Tf/dTout) ! number of output times

  ! parameters that can be modified via CLI args or are derived
  real(c_double) :: hs         = 0.01d0   ! slow step size
  real(c_double) :: e          = 0.5d0    ! fast/slow coupling strength
  real(c_double) :: G          = -100.0d0 ! stiffness at slow time scale
  real(c_double) :: w          = 100.0d0  ! time-scale separation factor
  real(c_double) :: reltol     = 0.01d0   ! integrator tolerances
  real(c_double) :: abstol     = 1e-11

contains

  !
  ! RHS functions
  !

  ! ff routine to compute the fast portion of the ODE RHS.
  integer(c_int) function ff(t, yvec, ydotvec, user_data) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: yvec
    type(N_Vector) :: ydotvec
    type(c_ptr)    :: user_data

    real(c_double) :: u, v
    real(c_double) :: tmp1, tmp2
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: ydotarr(:)

    ! get N_Vector data arrays
    yarr => FN_VGetArrayPointer(yvec)
    ydotarr => FN_VGetArrayPointer(ydotvec)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the RHS function:
    !  [0  0]*[(-1+u^2-r(t))/(2*u)] + [         0          ]
    !  [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))]
    tmp1 = (-ONE+u*u-r(t))/(TWO*u)
    tmp2 = (-TWO+v*v-s(t))/(TWO*v)
    ydotarr(1) = ZERO
    ydotarr(2) = e*tmp1 - tmp2 + sdot(t)/(TWO*vtrue(t))

    ! return success
    ierr = 0
    return

  end function

  ! fs routine to compute the slow portion of the ODE RHS.
  integer(c_int) function fs(t, yvec, ydotvec, user_data) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: yvec
    type(N_Vector) :: ydotvec
    type(c_ptr)    :: user_data

    real(c_double) :: u, v
    real(c_double) :: tmp1, tmp2
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: ydotarr(:)

    ! get N_Vector data arrays
    yarr => FN_VGetArrayPointer(yvec)
    ydotarr => FN_VGetArrayPointer(ydotvec)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the RHS function:
    !  [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
    !  [0 0] [(-2+v^2-s(t))/(2*v)]    [      0      ]
    tmp1 = (-ONE+u*u-r(t))/(TWO*u)
    tmp2 = (-TWO+v*v-s(t))/(TWO*v)
    ydotarr(1) = G*tmp1 + e*tmp2 + rdot(t)/(TWO*u)
    ydotarr(2) = ZERO

    ! return success
    ierr = 0
    return

  end function

  ! fse routine to compute the slow portion of the ODE RHS.
  integer(c_int) function fse(t, yvec, ydotvec, user_data) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: yvec
    type(N_Vector) :: ydotvec
    type(c_ptr)    :: user_data

    real(c_double) :: u, v
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: ydotarr(:)

    ! get N_Vector data arrays
    yarr => FN_VGetArrayPointer(yvec)
    ydotarr => FN_VGetArrayPointer(ydotvec)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the slow explicit RHS function:
    !  [rdot(t)/(2*u)]
    !  [      0      ]
    ydotarr(1) = rdot(t)/(TWO*u)
    ydotarr(2) = ZERO

    ! return success
    ierr = 0
    return

  end function

  ! fsi routine to compute the slow portion of the ODE RHS.(currently same as fse)
  integer(c_int) function fsi(t, yvec, ydotvec, user_data) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: yvec
    type(N_Vector) :: ydotvec
    type(c_ptr)    :: user_data

    real(c_double) :: u, v
    real(c_double) :: tmp1, tmp2
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: ydotarr(:)

    ! get N_Vector data arrays
    yarr => FN_VGetArrayPointer(yvec)
    ydotarr => FN_VGetArrayPointer(ydotvec)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the slow implicit RHS function:
    !  [G e]*[(-1+u^2-r(t))/(2*u))]
    !  [0 0] [(-2+v^2-s(t))/(2*v)]
    tmp1 = (-ONE+u*u-r(t))/(TWO*u)
    tmp2 = (-TWO+v*v-s(t))/(TWO*v)
    ydotarr(1) = G*tmp1 + e*tmp2
    ydotarr(2) = ZERO

    ! return success
    ierr = 0
    return

  end function

  integer(c_int) function fn(t, yvec, ydotvec, user_data) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: yvec
    type(N_Vector) :: ydotvec
    type(c_ptr)    :: user_data

    real(c_double) :: u, v
    real(c_double) :: tmp1, tmp2
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: ydotarr(:)

    ! get N_Vector data arrays
    yarr => FN_VGetArrayPointer(yvec)
    ydotarr => FN_VGetArrayPointer(ydotvec)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the RHS function:
    !  [G e]*[(-1+u^2-r(t))/(2*u))] + [rdot(t)/(2*u)]
    !  [e -1] [(-2+v^2-s(t))/(2*v)]   [sdot(t)/(2*vtrue(t))]
    tmp1 = (-ONE+u*u-r(t))/(TWO*u)
    tmp2 = (-TWO+v*v-s(t))/(TWO*v)
    ydotarr(1) = G*tmp1 + e*tmp2 + rdot(t)/(TWO*u)
    ydotarr(2) = e*tmp1 - tmp2 + sdot(t)/(TWO*vtrue(t))

    ! return success
    ierr = 0
    return

  end function

  integer(c_int) function f0(t, yvec, ydotvec, user_data) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: yvec
    type(N_Vector) :: ydotvec
    type(c_ptr)    :: user_data

    call FN_VConst(ZERO, ydotvec)

    ! return success
    ierr = 0
    return

  end function

  !
  ! Jacobian functions
  !

  integer(c_int) function Js(t, y, fy, J, user_data, tmp1, tmp2, tmp3) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value  :: t
    type(N_Vector)  :: y
    type(N_Vector)  :: fy
    type(SUNMatrix) :: J
    type(c_ptr)     :: user_data
    type(N_Vector)  :: tmp1, tmp2, tmp3

    real(c_double) :: u, v
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: Jarr(:)

    ! get N_Vector data arrays
    yarr => FN_VGetArrayPointer(y)

    ! get Jacobian data array
    Jarr => FSUNDenseMatrix_Data(J)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the Jacobian:
    !  [G/2 + (w*(1+r(t))-rdot(t))/(2*u^2)   e/2 + e*(2+s(t))/(2*v^2)]
    !  [                 0                              0            ]
    Jarr(1) = G/TWO + (G*(ONE+r(t))-rdot(t))/(2*u*u)
    Jarr(2) = ZERO
    Jarr(3) = e/TWO + e*(TWO+s(t))/(TWO*v*v)
    Jarr(4) = ZERO

    ! return success
    ierr = 0
    return

  end function

  integer(c_int) function Jsi(t, y, fy, J, user_data, tmp1, tmp2, tmp3) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value :: t
    type(N_Vector)  :: y
    type(N_Vector)  :: fy
    type(SUNMatrix) :: J
    type(c_ptr)     :: user_data
    type(N_Vector)  :: tmp1, tmp2, tmp3

    real(c_double) :: u, v
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: Jarr(:)

    ! get N_Vector data array
    yarr => FN_VGetArrayPointer(y)

    ! get Jacobian data array
    Jarr => FSUNDenseMatrix_Data(J)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the Jacobian:
    !  [G/2 + (G*(1+r(t)))/(2*u^2)   e/2+e*(2+s(t))/(2*v^2)]
    !  [                 0                             0   ]
    Jarr(1) = G/TWO + (G*(ONE+r(t)))/(2*u*u)
    Jarr(2) = ZERO
    Jarr(3) = e/TWO + e*(TWO+s(t))/(TWO*v*v)
    Jarr(4) = ZERO

    ! return success
    ierr = 0
    return

  end function

  integer(c_int) function Jn(t, y, fy, J, user_data, tmp1, tmp2, tmp3) &
    result(ierr) bind(C)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value  :: t
    type(N_Vector)  :: y
    type(N_Vector)  :: fy
    type(SUNMatrix) :: J
    type(c_ptr)     :: user_data
    type(N_Vector)  :: tmp1, tmp2, tmp3

    real(c_double) :: u, v
    real(c_double), pointer :: yarr(:)
    real(c_double), pointer :: Jarr(:)

    ! get N_Vector data array
    yarr => FN_VGetArrayPointer(y)

    ! get Jacobian data array
    Jarr => FSUNDenseMatrix_Data(J)

    ! extract variables
    u = yarr(1)
    v = yarr(2)

    ! fill in the Jacobian:
    ! [G/2 + (G*(1+r(t))-rdot(t))/(2*u^2)     e/2 + e*(2+s(t))/(2*v^2)]
    ! [e/2 + e*(1+r(t))/(2*u^2)              -1/2 - (2+s(t))/(2*v^2)]
    Jarr(1) = G/TWO + (G*(1+r(t))-rdot(t))/(TWO*u*u)
    Jarr(2) = e/TWO + e*(ONE+r(t))/(TWO*u*u)
    Jarr(3) = e/TWO + e*(TWO+s(t))/(TWO*v*v)
    Jarr(4) = -ONE/TWO - (TWO+s(t))/(TWO*v*v)

    ! return success
    ierr = 0
    return

  end function

  ! Helper functions
  real(c_double) function r(t) &
    result(result)

    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double) :: t

    result = 0.5d0*cos(t)
    return

  end function

  real(c_double) function s(t) &
    result(result)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double) :: t

    result = cos(w*t)
    return

  end function

  real(c_double) function rdot(t) &
    result(result)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double) :: t

    result = -0.5d0*sin(t)
    return

  end function

  real(c_double) function sdot(t) &
    result(result)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double) :: t

    result = -w*sin(w*t)
    return

  end function

  real(c_double) function utrue(t) &
    result(result)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double) :: t

    result = sqrt(ONE+r(t))
    return

  end function

  real(c_double) function vtrue(t) &
    result(result)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double) :: t

    result = sqrt(TWO+s(t))
    return

  end function

  integer(c_int) function Ytrue(t, y) &
    result(ierr)

    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double), value :: t
    type(N_Vector) :: y

    real(c_double), pointer :: yarr(:)

    yarr => FN_VGetArrayPointer(y)

    yarr(1) = utrue(t)
    yarr(2) = vtrue(t)

    ierr = 0
    return

  end function

end module

program main

  use ode_mod
  implicit none


  ! general problem variables
  type(c_ptr) sunctx                       ! SUNDIALS simulation context
  integer(c_int) :: retval                 ! reusable error-checking flag
  type(N_Vector), pointer :: y             ! vector for the solution
  real(c_double), pointer :: yarr(:)       ! array of data for y vector
  type(c_ptr) :: arkode_mem                ! ARKODE memory structure
  type(c_ptr) :: inner_arkode_mem          ! ARKODE memory structure
  type(c_ptr) :: inner_stepper             ! inner stepper
  type(c_ptr) :: BTF, BTs                  ! fast/slow method Butcher table
  type(c_ptr) :: SC                        ! slow coupling coefficients
  type(SUNMatrix), pointer :: MATf         ! matrix for fast solver
  type(SUNLinearSolver), pointer :: LSf    ! fast linear solver object
  type(SUNMatrix), pointer :: MATs         ! matrix for slow solver
  type(SUNLinearSolver), pointer :: LSs    ! slow linear solver object
  integer(c_int) :: solve_type = 0         ! problem configuration type
  logical :: implicit_slow
  logical :: imex_slow = .FALSE.
  real(c_double) :: hf, gamma, beta, t, tret(1), tout
  real(c_double) :: uerr, verr, uerrtot, verrtot, errtot
  real(c_double), allocatable :: Af(:,:), bf(:), cf(:), df(:) ! Arrays for fast Butcher table, NOTE: must be in row-major order
  real(c_double), allocatable :: As(:,:), bs(:), cs(:), ds(:) ! Arrays for slow Butcher table, NOTE: must be in row-major order
  integer(c_int) :: iout, argc, argi
  integer(c_long) :: nsts(1), nstf(1), nfse(1), nfsi(1), nff(1)
  integer(c_long) :: nnif(1), nncf(1), njef(1), nnis(1), nncs(1), njes(1), tmp(1)
  character(len=32), dimension(:), allocatable :: argv

  !
  ! Initialization
  !

  arkode_mem = c_null_ptr
  inner_arkode_mem = c_null_ptr
  inner_stepper = c_null_ptr
  BTs = c_null_ptr
  SC = c_null_ptr
  MATf => null()
  MATs => null()
  LSf => null()
  LSs => null()

  argc = command_argument_count()
  allocate(argv(argc))  ! I've omitted checking the return status of the allocation

  do argi = 1, argc
    call get_command_argument(argi, argv(argi))
  end do

  ! Retrieve the command-line options: solve_type h G w e */
  if (argc > 0) read(argv(1), *) solve_type
  if (argc > 1) read(argv(2), *) hs
  if (argc > 2) read(argv(3), *) G
  if (argc > 3) read(argv(4), *) w
  if (argc > 4) read(argv(5), *) e

  ! Check arguments for validity
  !   0 <= solve_type <= 9
  !   G < 0.0
  !   h > 0
  !   h < 1/|G| (explicit slow)
  !   w >= 1.0
  if ((solve_type < 0) .or. (solve_type > 9)) then
    print *, "ERROR: solve_type be an integer in [0,9]"
    stop -1
  end if
  if (G >= ZERO) then
    print *, "ERROR: G must be a negative real number"
    stop -1
  end if
  implicit_slow = .false.
  if ((solve_type == 4) .or. (solve_type == 7)) then
    implicit_slow = .true.
  end if
  if ((solve_type == 8) .or. (solve_type == 9)) then
    implicit_slow = .true.
    imex_slow = .true.
  end if
  if (hs <= ZERO) then
    print *, "ERROR: hs must be in positive"
    stop -1
  end if
  if ((hs > ONE/abs(G)) .and. (.not. implicit_slow)) then
    print *, "ERROR: hs must be in (0, 1/|G|)"
    stop -1
  end if
  if (w < ONE) then
    print *, "ERROR: w must be >= 1.0"
    stop -1
  end if
  hf = hs/w

  ! Initial problem output (and set implicit solver tolerances as needed)
  print *, "Multirate nonlinear Kvaerno-Prothero-Robinson test problem:"
  print '(A,F4.2,A,F4.2,A)', "    time domain: (", T0, ", ", Tf, "]"
  print '(A,E12.4,A)', "    hs = ", hs
  print '(A,E12.4,A)', "    hf = ", hf
  print '(A,E12.4,A)', "     G = ", G
  print '(A,E12.4,A)', "     w = ", w
  print '(A,E12.4,A)', "     e = ", e

  select case (solve_type)
    case (0)
      print *, "   solver: exp-3/exp-3 (standard MIS)"
    case (1)
      print *, "   solver: none/exp-3 (no slow, explicit fast)"
    case (2)
      reltol = max(hs*hs*hs, real(1e-10,8))
      abstol = 1e-11
      print *, "   solver: none/dirk-3 (no slow, dirk fast)"
      print '(A,E12.4,A,E12.4)', "    reltol: ", reltol, " abstol: ", abstol
    case (3)
      print *, "   solver: exp-3/none (explicit slow, no fast)"
    case (4)
      reltol = max(hs*hs, real(1e-10,8))
      abstol = 1e-11
      print *, "   solver: dirk-2/none (dirk slow, no fast)"
      print '(A,E12.4,A,E12.4)', "    reltol: ", reltol, " abstol: ", abstol
    case (5)
      print *, "   solver: exp-4/exp-4 (MRI-GARK-ERK45a / ERK-4-4)"
    case (6)
      print *, "   solver: exp-4/exp-3 (MRI-GARK-ERK45a / ERK-3-3)"
    case (7)
      reltol = max(hs*hs*hs, real(1e-10,8))
      abstol = 1e-11
      print *, "   solver: dirk-3/exp-3 (MRI-GARK-ESDIRK34a / ERK-3-3) -- solve decoupled"
      print '(A,E12.4,A,E12.4)', "    reltol: ", reltol, " abstol: ", abstol
    case (8)
      reltol = max(hs*hs*hs, real(1e-10,8))
      abstol = 1e-11
      print *, "   solver: ars343/exp-3 (IMEX-MRI3b / ERK-3-3) -- solve decoupled"
      print '(A,E12.4,A,E12.4)', "    reltol: ", reltol, " abstol: ", abstol
    case (9)
      reltol = max(hs*hs*hs*hs, real(1e-14,8))
      abstol = 1e-14
      print *, "   solver: imexark4/exp-4 (IMEX-MRI4 / ERK-4-4) -- solve decoupled"
      print '(A,E12.4,A,E12.4)', "    reltol: ", reltol, " abstol: ", abstol
  end select

  ! Create the SUNDIALS context object for this simulation
  retval = FSUNContext_Create(c_null_ptr, sunctx)
  call check_retval(retval, 'FSUNContext_Create')

  ! Create and initialize serial vector for the solution
  y => FN_VNew_Serial(NEQ, sunctx)
  if (.not. associated(y)) then
    print *, 'ERROR: N_VNew_Serial failed'
    stop 1
  end if
  yarr => FN_VGetArrayPointer(y)

  retval = Ytrue(T0, y)
  call check_retval(retval, 'Ytrue')

  !
  ! Create the fast integrator and set options
  !

  ! Initialize the fast integrator. Specify the fast right-hand side
  ! function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the inital time T0,
  ! and the initial dependent variable vector y.

  if (solve_type == 0 .or. solve_type == 6 .or. solve_type == 7 .or. solve_type == 8) then
    ! erk-3-3 fast solver
    inner_arkode_mem = FARKStepCreate(c_funloc(ff), c_null_funptr, T0, y, sunctx)
    allocate(Af(3,3))
    allocate(bf(3))
    allocate(cf(3))
    allocate(df(3))
    Af = 0.d0
    bf = 0.d0
    cf = 0.d0
    df = 0.d0
    Af(1,2) = 0.5d0
    Af(1,3) = -ONE
    Af(2,3) = TWO
    bf(1) = ONE/6.0d0
    bf(2) = TWO/3.0d0
    bf(3) = ONE/6.0d0
    df(3) = ONE
    cf(2) = 0.5d0
    cf(3) = ONE
    BTf = FARKodeButcherTable_Create(3, 3, 2, cf, Af, bf, df)
    retval = FARKStepSetTables(inner_arkode_mem, 3, 2, c_null_ptr, BTf)
    call check_retval(retval, "FARKStepSetTables")
  else if (solve_type == 1) then
    ! erk-3-3 fast solver (full problem)
    inner_arkode_mem = FARKStepCreate(c_funloc(fn), c_null_funptr, T0, y, sunctx)
    allocate(Af(3,3))
    allocate(bf(3))
    allocate(cf(3))
    allocate(df(3))
    Af = 0.d0
    bf = 0.d0
    cf = 0.d0
    df = 0.d0
    Af(1,2) = 0.5d0
    Af(1,3) = -ONE
    Af(2,3) = TWO
    bf(1) = ONE/6.0d0
    bf(2) = TWO/3.0d0
    bf(3) = ONE/6.0d0
    df(3) = ONE
    cf(2) = 0.5d0
    cf(3) = ONE
    BTf = FARKodeButcherTable_Create(3, 3, 2, cf, Af, bf, df)
    retval = FARKStepSetTables(inner_arkode_mem, 3, 2, c_null_ptr, BTf)
    call check_retval(retval, "FARKStepSetTables")
  else if (solve_type == 5 .or. solve_type == 9) then
    ! erk-4-4 fast solver
    inner_arkode_mem = FARKStepCreate(c_funloc(ff), c_null_funptr, T0, y, sunctx)
    allocate(Af(4,4))
    allocate(bf(4))
    allocate(cf(4))
    allocate(df(4))
    Af = 0.d0
    bf = 0.d0
    cf = 0.d0
    df = 0.d0
    Af(1,2) = 0.5d0
    Af(2,3) = 0.5d0
    Af(3,4) = ONE
    bf(1) = ONE/6.0d0
    bf(2) = ONE/3.0d0
    bf(3) = ONE/3.0d0
    bf(4) = ONE/6.0d0
    cf(2) = 0.5d0
    cf(3) = 0.5d0
    cf(4) = ONE
    BTf = FARKodeButcherTable_Create(4, 4, 0, cf, Af, bf, df)
    retval = FARKStepSetTables(inner_arkode_mem, 4, 0, c_null_ptr, BTf)
    call check_retval(retval, "FARKStepSetTables")
  else if (solve_type == 2) then
    ! esdirk-3-3 fast solver (full problem)
    inner_arkode_mem = FARKStepCreate(c_null_ptr, c_funloc(fn), T0, y, sunctx)
    beta  = sqrt(3.0d0)/6.0d0 + 0.5d00
    gamma = (-ONE/8.0d0)*(sqrt(3.0d0)+ONE)
    allocate(Af(3,3))
    allocate(bf(3))
    allocate(cf(3))
    allocate(df(3))
    Af = 0.d0
    bf = 0.d0
    cf = 0.d0
    df = 0.d0
    Af(1,2) = 4.0d0*gamma+TWO*beta
    Af(2,2) = ONE-4.0d0*gamma-TWO*beta
    Af(1,3) = 0.5d0-beta-gamma
    Af(2,3) = gamma
    Af(3,3) = beta
    bf(1) = ONE/6.0d0
    bf(2) = ONE/6.0d0
    bf(3) = TWO/3.0d0
    cf(2) = ONE
    cf(3) = 0.5d0
    BTf = FARKodeButcherTable_Create(3, 3, 0, cf, Af, bf, df)
    retval = FARKStepSetTables(inner_arkode_mem, 3, 0, BTf, c_null_ptr)
    call check_retval(retval, "FARKStepSetTables")
    MATf => FSUNDenseMatrix(NEQ, NEQ, sunctx)
    LSf => FSUNLinSol_Dense(y, MATf, sunctx)
    retval = FARKStepSetLinearSolver(inner_arkode_mem, LSf, MATf)
    call check_retval(retval, "FARKStepSetLinearSolver")
    retval = FARKStepSetJacFn(inner_arkode_mem, c_funloc(Jn))
    call check_retval(retval, "FARKStepSetJacFn")
    retval = FARKStepSStolerances(inner_arkode_mem, reltol, abstol)
    call check_retval(retval, "FARKStepSStolerances")
  else if (solve_type == 3 .or. solve_type == 4) then
    ! no fast dynamics ('evolve' explicitly w/ erk-3-3)
    inner_arkode_mem = FARKStepCreate(c_funloc(f0), c_null_funptr, T0, y, sunctx)
    allocate(Af(3,3))
    allocate(bf(3))
    allocate(cf(3))
    allocate(df(3))
    Af = 0.d0
    bf = 0.d0
    cf = 0.d0
    df = 0.d0
    Af(1,2) = 0.5d0
    Af(1,3) = -ONE
    Af(2,3) = TWO
    bf(1) = ONE/6.0d0
    bf(2) = TWO/3.0d0
    bf(3) = ONE/6.0d0
    df(1) = ONE
    cf(2) = 0.5d0
    cf(3) = ONE
    BTf = FARKodeButcherTable_Create(3, 3, 2, cf, Af, bf, df)
    retval = FARKStepSetTables(inner_arkode_mem, 3, 2, c_null_ptr, BTf)
    call check_retval(retval, "FARKStepSetTables")
  end if

  if (.not. c_associated(inner_arkode_mem)) then
    print *, 'ERROR: inner_arkode_mem = NULL'
    stop 1
  end if

  ! Set the fast step size */
  retval = FARKStepSetFixedStep(inner_arkode_mem, hf)
  if (retval /= 0) then
    print *, 'ERROR: FARKStepSetFixedStep failed'
    stop 1
  end if

  ! Create inner stepper */
  retval = FARKStepCreateMRIStepInnerStepper(inner_arkode_mem, inner_stepper)
  if (retval /= 0) then
    print *, 'ERROR: FARKStepCreateMRIStepInnerStepper failed'
    stop 1
  end if

  !
  ! Create the slow integrator and set options
  !

  ! Initialize the slow integrator. Specify the slow right-hand side
  ! function in y'=fs(t,y)+ff(t,y) = fse(t,y)+fsi(t,y)+ff(t,y), the inital time
  ! T0, the initial dependent variable vector y, and the fast integrator.
  if (solve_type == 0) then
    ! KW3 slow solver
    arkode_mem = FMRIStepCreate(c_funloc(fs), c_null_funptr, T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_MIS_KW3)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
  else if (solve_type == 3) then
    ! KW3 slow solver (full problem)
    arkode_mem = FMRIStepCreate(c_funloc(fn), c_null_funptr, T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_MIS_KW3)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
  else if (solve_type == 5 .or. solve_type == 6) then
    ! MRI-GARK-ERK45a slow solver
    arkode_mem = FMRIStepCreate(c_funloc(fs), c_null_funptr, T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ERK45a)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
  else if (solve_type == 1 .or. solve_type == 2) then
    ! no slow dynamics (use ERK-2-2)
    arkode_mem = FMRIStepCreate(c_funloc(f0), c_null_funptr, T0, y, inner_stepper, sunctx)
    if (.not. c_associated(arkode_mem)) then
      print *, 'ERROR: arkode_mem = NULL'
      stop 1
    end if
    allocate(As(2,2))
    allocate(bs(2))
    allocate(cs(2))
    allocate(ds(2))
    As = 0.d0
    bs = 0.d0
    cs = 0.d0
    ds = 0.d0
    As(1,2) = TWO/3.0d0
    bs(1) = 0.25d0
    bs(2) = 0.75d0
    cs(2) = TWO/3.0d0
    BTs = FARKodeButcherTable_Create(2, 2, 0, cs, As, bs, ds)
    SC = FMRIStepCoupling_MIStoMRI(BTs, 2, 0)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
  else if (solve_type == 4) then
    ! dirk-2 (trapezoidal), solve-decoupled slow solver
    arkode_mem = FMRIStepCreate(c_null_ptr, c_funloc(fn), T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_MRI_GARK_IRK21a)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
    MATs => FSUNDenseMatrix(NEQ, NEQ, sunctx)
    LSs => FSUNLinSol_Dense(y, MATs, sunctx)
    retval = FMRIStepSetLinearSolver(arkode_mem, LSs, MATs)
    call check_retval(retval, "FMRIStepSetLinearSolver")
    retval = FMRIStepSetJacFn(arkode_mem, c_funloc(Jn))
    call check_retval(retval, "FMRIStepSetJacFn")
    retval = FMRIStepSStolerances(arkode_mem, reltol, abstol)
    call check_retval(retval, "FMRIStepSStolerances")
  else if (solve_type == 7) then
    ! MRI-GARK-ESDIRK34a, solve-decoupled slow solver
    arkode_mem = FMRIStepCreate(c_null_ptr, c_funloc(fs), T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_MRI_GARK_ESDIRK34a)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
    MATs => FSUNDenseMatrix(NEQ, NEQ, sunctx)
    LSs => FSUNLinSol_Dense(y, MATs, sunctx)
    retval = FMRIStepSetLinearSolver(arkode_mem, LSs, MATs)
    call check_retval(retval, "FMRIStepSetLinearSolver")
    retval = FMRIStepSetJacFn(arkode_mem, c_funloc(Js))
    call check_retval(retval, "FMRIStepSetJacFn")
    retval = FMRIStepSStolerances(arkode_mem, reltol, abstol)
    call check_retval(retval, "FMRIStepSStolerances")
  else if (solve_type == 8) then
    ! IMEX-MRI-GARK3b, solve-decoupled slow solver
    arkode_mem = FMRIStepCreate(c_funloc(fse), c_funloc(fsi), T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK3b)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    call check_retval(retval, "FMRIStepSetCoupling")
    MATs => FSUNDenseMatrix(NEQ, NEQ, sunctx)
    LSs => FSUNLinSol_Dense(y, MATs, sunctx)
    retval = FMRIStepSetLinearSolver(arkode_mem, LSs, MATs)
    call check_retval(retval, "FMRIStepSetLinearSolver")
    retval = FMRIStepSetJacFn(arkode_mem, c_funloc(Jsi))
    call check_retval(retval, "FMRIStepSetJacFn")
    retval = FMRIStepSStolerances(arkode_mem, reltol, abstol)
    call check_retval(retval, "FMRIStepSStolerances")
  else if (solve_type == 9) then
    ! IMEX-MRI-GARK4, solve-decoupled slow solver
    arkode_mem = FMRIStepCreate(c_funloc(fse), c_funloc(fsi), T0, y, inner_stepper, sunctx)
    SC = FMRIStepCoupling_LoadTable(ARKODE_IMEX_MRI_GARK4)
    retval = FMRIStepSetCoupling(arkode_mem, SC)
    MATs => FSUNDenseMatrix(NEQ, NEQ, sunctx)
    LSs => FSUNLinSol_Dense(y, MATs, sunctx)
    retval = FMRIStepSetLinearSolver(arkode_mem, LSs, MATs)
    call check_retval(retval, "FMRIStepSetLinearSolver")
    retval = FMRIStepSetJacFn(arkode_mem, c_funloc(Jsi))
    call check_retval(retval, "FMRIStepSetJacFn")
    retval = FMRIStepSStolerances(arkode_mem, reltol, abstol)
    call check_retval(retval, "FMRIStepSStolerances")
  end if

  if (.not. c_associated(arkode_mem)) then
    print *, 'ERROR: arkode_mem = NULL'
    stop 1
  end if

  ! Set the slow step size
  retval = FMRIStepSetFixedStep(arkode_mem, hs)
  call check_retval(retval, "FMRIStepSetFixedStep")

  !
  ! Integrate ODE
  !

  ! Main time-stepping loop: calls MRIStepEvolve to perform the
  ! integration, then prints results. Stops when the final time
  ! has been reached
  t = T0
  tout = T0+dTout
  uerr = ZERO
  verr = ZERO
  uerrtot = ZERO
  verrtot = ZERO
  errtot  = ZERO
  print *, "        t           u           v       uerr      verr"
  print *, "   ------------------------------------------------------"
  print '(A, F10.6, A, F10.6, A, F10.6, A, E8.2, A, E8.2)', &
    "   ", t, "  ", yarr(1), "  ", yarr(2), "  ", uerr, "  ", verr

  do iout = 1, Nt
    ! call integrator
    retval = FMRIStepEvolve(arkode_mem, tout, y, tret, ARK_NORMAL)
    call check_retval(retval, "FMRIStepEvolve")

    ! access/print solution and error
    uerr = abs(yarr(1)-utrue(tret(1)))
    verr = abs(yarr(2)-vtrue(tret(1)))
    print '(A, F10.6, A, F10.6, A, F10.6, A, E8.2, A, E8.2)', &
    "   ", tret(1), "  ", yarr(1), "  ", yarr(2), "  ", uerr, "  ", verr
    uerrtot = uerrtot + uerr*uerr
    verrtot = verrtot + verr*verr
    errtot = errtot + uerr*uerr + verr*verr

    ! successful solve: update time
    tout = tout + dTout
    if (tout > Tf) then
      tout = Tf
    end if
  enddo
  uerrtot = sqrt(uerrtot / Nt)
  verrtot = sqrt(verrtot / Nt)
  errtot = sqrt(errtot / Nt / 2)
  print *, "   ------------------------------------------------------"

  !
  ! Finalize
  !

  ! Get some slow integrator statistics
  retval = FMRIStepGetNumSteps(arkode_mem, nsts)
  retval = FMRIStepGetNumRhsEvals(arkode_mem, nfse, nfsi)

  ! Get some fast integrator statistics
  retval = FARKStepGetNumSteps(inner_arkode_mem, nstf)
  retval = FARKStepGetNumRhsEvals(inner_arkode_mem, nff, tmp)

  ! Print some final statistics
  print *, "Final Solver Statistics:"
  print '(A, I7, A, I7)', "   Steps: nsts = ", nsts, ", nstf = ", nstf
  print '(A, E9.3, A, E9.3, A, E9.3)', "   u error = ", uerrtot, ", v error = ", verrtot, ", total error = ", errtot
  if (imex_slow) then
    print '(A, I7, A, I7, A, I7)', "   Total RHS evals: Fse = ", nfse(1), ", Fsi = ", nfsi(1), ", Ff = ", nff(1)
  else if (implicit_slow) then
    print '(A, I7, A, I7)', "   Total RHS evals: Fs = ", nfsi(1), ", Ff = ", nff(1)
  else
    print '(A, I7, A, I7)', "   Total RHS evals: Fs = ", nfse(1), ", Ff = ", nff(1)
  end if

  ! Get/print slow integrator decoupled implicit solver statistics
  if ((solve_type == 4) .or. (solve_type == 7) .or. &
      (solve_type == 8) .or. (solve_type == 9)) then
    retval = FMRIStepGetNonlinSolvStats(arkode_mem, nnis, nncs)
    call check_retval(retval, "MRIStepGetNonlinSolvStats")
    retval = FMRIStepGetNumJacEvals(arkode_mem, njes)
    call check_retval(retval, "MRIStepGetNumJacEvals")
    print '(A, I7)', "   Slow Newton iters = ", nnis
    print '(A, I7)', "   Slow Newton conv fails = ", nncs
    print '(A, I7)', "   Slow Jacobian evals = ", njes
  end if

  ! Get/print fast integrator implicit solver statistics
  if (solve_type == 2) then
    retval = FARKStepGetNonlinSolvStats(inner_arkode_mem, nnif, nncf)
    call check_retval(retval, "ARKStepGetNonlinSolvStats")
    retval = FARKStepGetNumJacEvals(inner_arkode_mem, njef)
    call check_retval(retval, "ARKStepGetNumJacEvals")
    print '(A, I7)', "   Fast Newton iters = ", nnif
    print '(A, I7)', "   Fast Newton conv fails = ", nncf
    print '(A, I7)', "   Fast Jacobian evals = ", njef
  end if

  ! Clean up and return
  if (allocated(Af)) deallocate(Af)
  if (allocated(bf)) deallocate(bf)
  if (allocated(cf)) deallocate(cf)
  if (allocated(df)) deallocate(df)
  if (allocated(As)) deallocate(As)
  if (allocated(bs)) deallocate(bs)
  if (allocated(cs)) deallocate(cs)
  if (allocated(ds)) deallocate(ds)
  call FN_VDestroy(y)                                ! Free y vector
  call FARKodeButcherTable_Free(BTs)                 ! Butcher table
  call FMRIStepCoupling_Free(SC)                     ! Free coupling coefficients
  if (associated(MATf)) call FSUNMatDestroy(MATf)    ! Free fast matrix
  if (associated(LSf)) retval = FSUNLinSolFree(LSf)  ! Free fast linear solver
  if (associated(MATs)) call FSUNMatDestroy(MATs)    ! Free slow matrix
  if (associated(LSs)) retval = FSUNLinSolFree(LSs)  ! Free slow linear solver
  call FARKStepFree(inner_arkode_mem)                ! Free fast integrator memory
  retval = FMRIStepInnerStepper_Free(inner_stepper)  ! Free inner stepper
  call FMRIStepFree(arkode_mem)                      ! Free slow integrator memory
  retval = FSUNContext_Free(sunctx)                  ! Free context

end program

subroutine check_retval(retval, name)
  use, intrinsic :: iso_c_binding

  character(len=*) :: name
  integer(c_int)   :: retval

  if (retval /= 0) then
    write(*,'(A,A,A)') 'ERROR: ', name,' returned nonzero'
    stop 1
  end if
end subroutine
