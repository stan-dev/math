! -----------------------------------------------------------------------------
! Programmer(s): David J. Gardner, Cody J. Balos @ LLNL
! -----------------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------------------
! This demonstration problem simulates the advection and reaction of three
! chemical species, u, v, and w, in a one dimensional domain. The reaction
! mechanism is a variation of the Brusselator problem from chemical kinetics.
! This is a PDE system with 3 components, Y = [u,v,w], satisfying the
! equations,
!
!    u_t = -c * u_x + A - (w+1) * u + v * u^2
!    v_t = -c * v_x + w * u - v * u^2
!    w_t = -c * w_x + (B - w) / ep - w * u
!
! for t in [0, 10], x in [0, xmax] with periodic boundary conditions. The
! initial condition is a Gaussian pertubation of the steady state
! solution without advection
!
!    u(0,x) = k1 * A / k4 + p(x)
!    v(0,x) = k2 * k4 * B / (k1 * k3 * A) + p(x)
!    w(0,x) = 3.0 + p(x)
!    p(x)   = alpha * e^( -(x - mu)^2 / (2*sigma^2) ).
!
! where alpha = 0.1, mu = xmax / 2.0, and sigma = xmax / 4.0.
! The reaction rates are set so k_1 = k_2 = k_3 = k_4 = k, and k_5 = k_6 =
! 1/5e-6. The spatial derivatives are discretized with first-order upwind
! finite differences. An IMEX method is used to evolve the system in time with
! the advection terms treated explicitly and the reaction terms implicitly. As
! the reactions are purely local, the code uses a custom nonlinear solver to
! exploit this locality so no parallel communication is needed in the implicit
! solves. NOUT outputs are printed at equal intervals, and run statistics are
! printed at the end.
!
! Command line options:
!  --help           prints help message
!  --printtime      print timing information
!  --monitor        print solution information to screen (slower)
!  --nout <int>     the number of output times
!  --nx <int>       number of spatial mesh intervals
!  --xmax <double>  maximum x value
!  --explicit       use explicit method instead of IMEX
!  --order <int>    method order
!  --global-nls     use a global newton nonlinear solver instead of task-local
!  --tf <double>    final time
!  --A <double>     A parameter value
!  --B <double>     B parameter value
!  --k <double>     reaction rate
!  --c <double>     advection speed
!  --rtol <double>  relative tolerance
!  --atol <double>  absolute tolerance
! --------------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod

  !======= Declarations =========
  implicit none
  save

  ! Number of chemical species
  integer, parameter :: Nvar = 3

  ! MPI variables
  integer :: comm   ! communicator
  integer :: myid   ! process ID
  integer :: nprocs ! total number of processes
  integer :: reqS   ! MPI send request handle
  integer :: reqR   ! MPI receive request handle

  ! Excahnge buffers
  real(c_double) :: Wsend(Nvar), Wrecv(Nvar)
  real(c_double) :: Esend(Nvar), Erecv(Nvar)

  ! Problem settings
  integer :: Nx   ! number of intervals (global)
  integer :: Npts ! number of spatial nodes (local)
  integer :: Neq  ! number of equations (local)

  double precision :: xmax ! maximum x value
  double precision :: dx   ! mesh spacing

  ! Problem parameters
  double precision :: c  ! advection speed
  double precision :: A  ! constant concentrations
  double precision :: B
  double precision :: k1 ! reaction rates
  double precision :: k2
  double precision :: k3
  double precision :: k4
  double precision :: k5
  double precision :: k6

  ! integrator options
  real(c_double) :: t0, tf     ! initial and final time
  real(c_double) :: rtol, atol ! relative and absolute tolerance
  integer(c_int) :: order      ! method order
  type(c_ptr)    :: DFID       ! ARKODE diagnostics file

  logical :: explicit  ! use explicit or IMEX method
  logical :: global    ! use global or task-local nonlinear solver
  logical :: fused     ! use fused vector operations
  logical :: monitor   ! enable diagnostic output
  logical :: printtime ! print timing

  ! output settings and variables
  integer                 :: nout    ! number of outputs
  type(N_Vector), pointer :: umask_s ! mask vectors for RMS computation
  type(N_Vector), pointer :: umask
  type(N_Vector), pointer :: vmask
  type(N_Vector), pointer :: wmask

contains

  ! --------------------------------------------------------------
  ! Right Hand Side (RHS) Functions
  ! --------------------------------------------------------------

  ! Compute the advection term
  integer(c_int) function Advection(t, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t          ! current time
    type(N_Vector)        :: sunvec_y   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    type(c_ptr)           :: user_data  ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: ydata(:)
    real(c_double), pointer :: fdata(:)

    ! local variables
    integer          :: i, j, idx1, idx2 ! loop counters and array indices
    double precision :: tmp              ! temporary value

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    ydata => FN_VGetArrayPointer(sunvec_y)
    fdata => FN_VGetArrayPointer(sunvec_f)

    ! Set output to zero
    fdata = 0.0d0

    ! Begin exchanging boundary information
    call ExchangeAllStart(sunvec_y)

    ! Iterate over domain interior, computing advection
    tmp = -c / dx

    if (c > 0.0d0) then

       ! right moving flow
       do j = 2,Npts
          do i = 1,Nvar
             idx1 = i + (j - 1) * Nvar
             idx2 = i + (j - 2) * Nvar
             fdata(idx1) = tmp * (ydata(idx1) - ydata(idx2))
          end do
       end do

    else if (c < 0.0d0) then

       ! left moving flow
       do j = 1,Npts - 1
          do i = 1,Nvar
             idx1 = i + (j - 1) * Nvar
             idx2 = i + j * Nvar
             fdata(idx1) = tmp * (ydata(idx2) - ydata(idx1))
          end do
       end do

    end if

    ! finish exchanging boundary information
    call ExchangeAllEnd()

    ! compute advection at local boundaries
    if (c > 0.0d0) then

       ! right moving flow (left boundary)
       fdata(1:Nvar) = tmp * (ydata(1:Nvar) - Wrecv)

    else if (c < 0.0) then

       ! left moving flow (right boundary)
       fdata(Nvar * Npts - 2 : Nvar * Npts) = &
            tmp * (Erecv - ydata(Nvar * Npts-2 : Nvar * Npts))

    end if

    ! return success
    ierr = 0

  end function Advection


  ! Compute the reaction term
  integer(c_int) function Reaction(t, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t          ! current time
    type(N_Vector)        :: sunvec_y   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    type(c_ptr)           :: user_data  ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: ydata(:)
    real(c_double), pointer :: fdata(:)

    ! local variables
    double precision :: u, v, w ! chemcial species
    integer          :: j, idx  ! loop counter and array index

    !======= Internals ============

    ! Get data arrays from SUNDIALS vectors
    ydata => FN_VGetArrayPointer(sunvec_y)
    fdata => FN_VGetArrayPointer(sunvec_f)

    ! iterate over domain, computing reactions
    if (explicit) then

       ! when integrating explicitly, we add to ydot as we expect it
       ! to hold the advection term already
       do j = 1,Npts

          idx = (j - 1) * Nvar

          u = ydata(idx + 1)
          v = ydata(idx + 2)
          w = ydata(idx + 3)

          fdata(idx + 1) = fdata(idx + 1) + k1 * A - k2 * w * u + k3 * u * u * v - k4 * u
          fdata(idx + 2) = fdata(idx + 2) + k2 * w * u - k3 * u * u * v
          fdata(idx + 3) = fdata(idx + 3) - k2 * w * u + k5 * B - k6 * w

       end do

    else

       ! set output to zero
       fdata = 0.0d0

       do j = 1,Npts

          idx = (j - 1) * Nvar

          u = ydata(idx + 1)
          v = ydata(idx + 2)
          w = ydata(idx + 3)

          fdata(idx + 1) = k1 * A - k2 * w * u + k3 * u * u * v - k4 * u
          fdata(idx + 2) = k2 * w * u - k3 * u * u * v
          fdata(idx + 3) = -k2 * w * u + k5 * B - k6 * w

       end do

    end if

    ! return success
    ierr = 0

  end function Reaction


  ! Compute the RHS as Advection + Reaction
  integer(c_int) function AdvectionReaction(t, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t          ! current time
    type(N_Vector)        :: sunvec_y   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    type(c_ptr)           :: user_data  ! user-defined data

    !======= Internals ============

    ! NOTE: The order in which Advection and Reaction are called is critical
    ! here. Advection must be computed first.

    ierr = Advection(t, sunvec_y, sunvec_f, user_data)
    if (ierr /= 0) return

    ierr = Reaction(t, sunvec_y, sunvec_f, user_data)
    if (ierr /= 0) return

    ! return success
    ierr = 0

  end function AdvectionReaction

end module ode_mod


module prec_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_matrix_mod
  use fsundials_linearsolver_mod

  !======= Declarations =========
  implicit none
  save

  ! preconditioner data
  type(SUNLinearSolver), pointer :: sunls_P  ! linear solver
  type(SUNMatrix),       pointer :: sunmat_P ! matrix

contains

  ! --------------------------------------------------------------
  ! Preconditioner functions
  ! --------------------------------------------------------------


  ! Sets P = I - gamma * J
  integer(c_int) function PSetup(t, sunvec_y, sunvec_f, jok, jcurPtr, gamma, &
       user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsunmatrix_dense_mod
    use fsundials_linearsolver_mod
    use fsunlinsol_dense_mod

    use ode_mod, only : Nvar, Npts, Neq, k2, k3, k4, k6

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t          ! current time
    type(N_Vector)        :: sunvec_y   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    integer(c_int), value :: jok        ! flag to signal for Jacobian update
    integer(c_int)        :: jcurPtr    ! flag to singal Jacobian is current
    real(c_double), value :: gamma      ! current gamma value
    type(c_ptr)           :: user_data  ! user-defined data

    ! local variables
    real(c_double), pointer :: ydata(:) ! vector data
    real(c_double), pointer :: pdata(:) ! matrix data

    double precision :: u, v, w        ! chemical species
    integer          :: i, idx, offset ! loop counter, array index, col offset

    !======= Internals ============

    ! access solution data
    ydata => FN_VGetArrayPointer(sunvec_y)
    pdata => FSUNDenseMatrix_Data(sunmat_P)

    ! update Jacobian
    if (jok == 0) then

       ! zero the matrix
       ierr = FSUNMatZero(sunmat_P)
       if (ierr /= 0) then
          print *, "Error: FSUNMatZero returned ",ierr
          return
       end if

       ! setup the block diagonal preconditioner matrix
       do i = 1,Npts

          ! set nodal value shortcuts
          idx = (i - 1) * Nvar

          u = ydata(idx + 1)
          v = ydata(idx + 2)
          w = ydata(idx + 3)

          ! fill in Jacobian entries for this mesh node

          ! first column (derivative with respect to u)
          offset = (i - 1) * Nvar * Neq + (i - 1) * Nvar

          pdata(offset + 1) = -k2 * w + 2.0d0 * k3 * u * v - k4
          pdata(offset + 2) =  k2 * w - 2.0d0 * k3 * u * v
          pdata(offset + 3) = -k2 * w

          ! second column (derivative with respect to v)
          offset = offset + Nvar * Npts

          pdata(offset + 1) =  k3 * u * u
          pdata(offset + 2) = -k3 * u * u
          pdata(offset + 3) =  0.0d0

          ! thrid column (derivative with respect to v)
          offset = offset + Neq

          pdata(offset + 1) = -k2 * u
          pdata(offset + 2) =  k2 * u
          pdata(offset + 3) = -k2 * u - k6

       end do

       ierr = FSUNMatScaleAddI(-gamma, sunmat_P)
       if (ierr /= 0) then
          print *, "Error: FSUNMatScaleAddI returned ",ierr
          return
       end if

       ! setup the linear system Pz = r
       ierr = FSUNLinSolSetup(sunls_P, sunmat_P)
       if (ierr /= 0) then
          print *, "Error: FSUNLinSolSetup returned ",ierr
          return
       end if

       ! indicate that J is now current
       jcurPtr = 1

    else

       jcurPtr = 0

    end if

    ! return success
    ierr = 0

  end function PSetup


  ! Solves Pz = r
  integer(c_int) function PSolve(t, sunvec_y, sunvec_f, sunvec_r, sunvec_z, &
       gamma, delta, lr, user_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_mpiplusx_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t          ! current time
    type(N_Vector)        :: sunvec_y   ! solution N_Vector
    type(N_Vector)        :: sunvec_f   ! rhs N_Vector
    type(N_Vector)        :: sunvec_r   ! rhs N_Vector
    type(N_Vector)        :: sunvec_z   ! rhs N_Vector
    real(c_double), value :: gamma      ! current gamma value
    real(c_double), value :: delta      ! current gamma value
    integer(c_int), value :: lr         ! left or right preconditioning
    type(c_ptr)           :: user_data  ! user-defined data

    ! shortcuts
    type(N_Vector), pointer :: z_local ! z vector data
    type(N_Vector), pointer :: r_local ! r vector data

    !======= Internals ============

    z_local => FN_VGetLocalVector_MPIPlusX(sunvec_z)
    r_local => FN_VGetLocalVector_MPIPlusX(sunvec_r)

    ! solve the task-local linear system Pz = r
    ierr = FSUNLinSolSolve(sunls_P, sunmat_P, z_local, r_local, delta)
    if (ierr /= 0) then
       print *, "Error: FSUNLinSolSolver returned ",ierr
       return
    end if

    ! return success
    ierr = 0

  end function PSolve

end module prec_mod


module nls_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod
  use fsundials_matrix_mod
  use fsundials_linearsolver_mod
  use fsundials_nonlinearsolver_mod

  !======= Declarations =========
  implicit none
  save

  ! ARKODE memory
  type(c_ptr), pointer :: arkmem

  ! task local nonlinear solver
  type(SUNNonlinearSolver), pointer :: sunnls_LOC

  ! nonlinear residual vectors
  type(c_ptr) :: zpred_ptr ! stage prediction vector
  type(c_ptr) :: z_ptr     ! current stage vector
  type(c_ptr) :: Fi_ptr    ! RHS vector
  type(c_ptr) :: sdata_ptr ! residual data

  ! node local linear solver and data
  type(N_Vector),        pointer :: sunvec_bnode ! node lobal rhs/solution vec
  type(SUNMatrix),       pointer :: sunmat_Jnode ! node local Jacobian
  type(SUNLinearSolver), pointer :: sunls_Jnode  ! node local linear solver

  ! nonlinear solver counters
  integer :: nnlfi    ! nonlinear function evals
  integer :: ncnf_loc ! nonlinear convergence fails

contains

  ! --------------------------------------------------------------
  ! (Non)linear system functions
  ! --------------------------------------------------------------

  integer(c_int) function TaskLocalNlsResidual(sunvec_zcor, sunvec_F, arkode_mem) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use farkode_arkstep_mod
    use fsundials_nvector_mod

    use ode_mod, only : Neq, Reaction

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_zcor ! current correction to predicted state
    type(N_Vector)     :: sunvec_F    ! nonlinear residual
    type(c_ptr), value :: arkode_mem  ! ARKODE memory structure

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: ycor_data(:)
    real(c_double), pointer :: F_data(:)

    ! residual data
    type(c_ptr) :: user_data

    real(c_double) :: tcur(1) ! current time
    real(c_double) :: gam(1)  ! current gamma

    ! SUNDIALS resiudal vectors
    type(N_Vector), pointer :: sunvec_zpred ! predicted stage vector
    type(N_Vector), pointer :: sunvec_z     ! current stage vector
    type(N_Vector), pointer :: sunvec_Fi    ! RHS vector
    type(N_Vector), pointer :: sunvec_sdata ! residual data vector

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: zpred_data(:)
    real(c_double), pointer :: z_data(:)
    real(c_double), pointer :: Fi_data(:)
    real(c_double), pointer :: sdata_data(:)

    integer :: i ! loop counter

    !======= Internals ============

    ! get nonlinear residual data
    ierr = FARKStepGetNonlinearSystemData(arkmem, tcur, zpred_ptr, z_ptr, &
         Fi_ptr, gam, sdata_ptr, user_data)
    if (ierr /= 0) then
       print *, "Error: FARKStepGetNonlinearSystemData returned ",ierr
       return
    end if

    ! get vectors from pointers
    sunvec_zpred => FN_VGetVecAtIndexVectorArray(zpred_ptr, 0)
    sunvec_z     => FN_VGetVecAtIndexVectorArray(z_ptr, 0)
    sunvec_Fi    => FN_VGetVecAtIndexVectorArray(Fi_ptr, 0)
    sunvec_sdata => FN_VGetVecAtIndexVectorArray(sdata_ptr, 0)

    ! get vector data arrays
    ycor_data  => FN_VGetArrayPointer(sunvec_zcor)
    F_data     => FN_VGetArrayPointer(sunvec_F)
    zpred_data => FN_VGetArrayPointer(sunvec_zpred)
    z_data     => FN_VGetArrayPointer(sunvec_z)
    Fi_data    => FN_VGetArrayPointer(sunvec_Fi)
    sdata_data => FN_VGetArrayPointer(sunvec_sdata)

    ! update "z" value as stored predictor + current corrector
    do i = 1,Neq
       z_data(i) = zpred_data(i) + ycor_data(i)
    end do

    ! compute implicit RHS and save for later
    ierr = Reaction(tcur(1), sunvec_z, sunvec_Fi, c_null_ptr)

    ! count calls to Fi as part of the nonlinear residual
    nnlfi = nnlfi + 1

    ! check RHS return value
    if (ierr /= 0) then
       print *, "Error: Reaction returned ",ierr
       return
    end if

    ! compute the nonlinear resiudal
    do i = 1,Neq
       F_data(i) = ycor_data(i) - sdata_data(i) - gam(1) * Fi_data(i)
    end do

    ! return success
    ierr = 0

  end function TaskLocalNlsResidual


  integer(c_int) function TaskLocalLSolve(sunvec_delta, arkode_mem) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use farkode_arkstep_mod
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsunmatrix_dense_mod
    use fsundials_linearsolver_mod

    use ode_mod, only : Nvar, Npts, k2, k3, k4, k6

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)     :: sunvec_delta ! input linear system rhs, ouput solution
    type(c_ptr), value :: arkode_mem   ! ARKODE memory structure

    ! residual data
    type(c_ptr) :: user_data

    real(c_double) :: tcur(1) ! current time
    real(c_double) :: gam(1)  ! current gamma

    type(N_Vector), pointer :: sunvec_z ! vector for evaluating J

    ! SUNDIALS vector data arrays
    real(c_double), pointer :: b_data(:)
    real(c_double), pointer :: z_data(:)
    real(c_double), pointer :: J_data(:)
    real(c_double), pointer :: bnode_data(:)

    double precision :: u, v, w ! chemical species
    integer          :: i, idx  ! loop counter and array index

    !======= Internals ============

    ! get nonlinear residual data
    ierr = FARKStepGetNonlinearSystemData(arkmem, tcur, zpred_ptr, z_ptr, &
         Fi_ptr, gam, sdata_ptr, user_data)
    if (ierr /= 0) then
       print *, "Error: FARKStepGetNonlinearSystemData returned ",ierr
       return
    end if

    ! get vectors from pointers
    sunvec_z => FN_VGetVecAtIndexVectorArray(z_ptr, 0)

    ! get data arrays
    b_data => FN_VGetArrayPointer(sunvec_delta)
    z_data => FN_VGetArrayPointer(sunvec_z)
    J_data => FSUNDenseMatrix_Data(sunmat_Jnode)

    bnode_data => FN_VGetArrayPointer(sunvec_bnode)

    ! solve the linear system at each mesh node
    do i = 1,Npts

       ! set nodal value shortcuts
       idx = (i - 1) * Nvar

       u = z_data(idx + 1)
       v = z_data(idx + 2)
       w = z_data(idx + 3)

       ! fill in Jacobian entries for this mesh node

       ! first column (derivative with respect to u)
       J_data(1) = -k2 * w + 2.0d0 * k3 * u * v - k4
       J_data(2) =  k2 * w - 2.0d0 * k3 * u * v
       J_data(3) = -k2 * w

       ! second column (derivative with respect to v)
       J_data(4) =  k3 * u * u
       J_data(5) = -k3 * u * u
       J_data(6) =  0.0d0

       ! thrid column (derivative with respect to v)
       J_data(7) = -k2 * u
       J_data(8) =  k2 * u
       J_data(9) = -k2 * u - k6

       ! I - gamma*J
       ierr = FSUNMatScaleAddI(-gam(1), sunmat_Jnode)
       if (ierr /= 0) then
          print *, "Error: FSUNMatScaleAddI returned ",ierr
          return
       end if

       ! grab just the portion of the vector "b" for this mesh node
       bnode_data = b_data(idx + 1 : idx + 3)

       ! setup the linear system
       ierr = FSUNLinSolSetup(sunls_Jnode, sunmat_Jnode)
       if (ierr /= 0) then
          print *, "Error: FSUNLinSolSolSetup returned ",ierr
          return
       end if

       ! solve the linear system
       ierr = FSUNLinSolSolve(sunls_Jnode, sunmat_Jnode, sunvec_bnode, &
            sunvec_bnode, 0.0d0)
       if (ierr /= 0) then
          print *, "Error: FSUNLinSolSolve returned ",ierr
          return
       end if

       ! set just the portion of the vector "b" for this mesh node
       b_data(idx + 1 : idx + 3) = bnode_data

    end do

    ! return success
    ierr = 0

  end function TaskLocalLSolve


  integer(SUNNonlinearSolver_Type) function TaskLocalNewton_GetType(sunnls) &
       result(id) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nonlinearsolver_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls ! nonlinear solver

    !======= Internals ============

    id = SUNNONLINEARSOLVER_ROOTFIND

  end function TaskLocalNewton_GetType


  integer(c_int) function TaskLocalNewton_Initialize(sunnls) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nonlinearsolver_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls ! nonlinear solver

    !======= Internals ============

    ! override default system and lsolve functions with local versions
    ierr = FSUNNonlinSolSetSysFn(sunnls_LOC, c_funloc(TaskLocalNlsResidual))
    if (ierr /= 0) then
       print *, "Error: FSUNNonlinSolSetSysFn returned ",ierr
       return
    end if

    ierr = FSUNNonlinSolSetLSolveFn(sunnls_LOC, c_funloc(TaskLocalLSolve))
    if (ierr /= 0) then
       print *, "Error: FSUNNonlinSolSetLSolveFn returned ",ierr
       return
    end if

    ierr = FSUNNonlinSolInitialize(sunnls_LOC)
    if (ierr /= 0) then
       print *, "Error: FSUNNonlinSolSetLSolveFn returned ",ierr
       return
    end if

    ! return success
    ierr = 0

  end function TaskLocalNewton_Initialize


  integer(c_int) function TaskLocalNewton_Solve(sunnls, sunvec_y0, sunvec_ycor, &
       sunvec_w, tol, callLSetup, arkode_mem) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_nonlinearsolver_mod
    use fnvector_mpiplusx_mod

    use ode_mod, only : comm

    !======= Declarations =========
    implicit none

    ! With MPI-3 use mpi_f08 is preferred
    include "mpif.h"

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls      ! nonlinear solver
    type(N_Vector)           :: sunvec_y0   ! initial guess
    type(N_Vector)           :: sunvec_ycor ! correction
    type(N_Vector)           :: sunvec_w    ! norm weights
    real(c_double), value    :: tol         ! solve tolerance
    integer(c_int), value    :: callLSetup  ! linear solver setup flag
    type(c_ptr)              :: arkode_mem  ! integrator memory

    ! local variables
    type(N_Vector), pointer :: sunvec_y0loc   ! node-local initial guess vector
    type(N_Vector), pointer :: sunvec_ycorloc ! node-local correction vector
    type(N_Vector), pointer :: sunvec_wloc    ! node-local weight vector

    integer :: solve_status, recover, nonrecover ! solve status, return value
    integer :: mpi_ierr                          ! MPI error status

    !======= Internals ============

    ! get MPI task local vectors
    sunvec_y0loc   => FN_VGetLocalVector_MPIPlusX(sunvec_y0)
    sunvec_ycorloc => FN_VGetLocalVector_MPIPlusX(sunvec_ycor)
    sunvec_wloc    => FN_VGetLocalVector_MPIPlusX(sunvec_w)

    ! each tasks solves the local nonlinear system
    ierr = FSUNNonlinSolSolve(sunnls_LOC, sunvec_y0loc, sunvec_ycorloc, &
         sunvec_wloc, tol, callLSetup, arkode_mem)
    solve_status = ierr

    ! if any process had a nonrecoverable failure, return it
    call MPI_Allreduce(solve_status, nonrecover, 1, MPI_INTEGER, MPI_MIN, comm, &
         mpi_ierr)
    ierr = nonrecover
    if (ierr < 0) return

    ! check if any process has a recoverable convergence failure and return
    ! success (recover == 0) or a recoverable error code (recover > 0)
    call MPI_Allreduce(solve_status, recover, 1, MPI_INTEGER, MPI_MAX, comm, &
         mpi_ierr)
    ierr = recover
    if (ierr /= 0) ncnf_loc = ncnf_loc + 1

  end function TaskLocalNewton_Solve


  integer(c_int) function TaskLocalNewton_Free(sunnls) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fsundials_linearsolver_mod
    use fsundials_nonlinearsolver_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls ! nonlinear solver

    !======= Internals ============

    ! free task-local solve structures
    call FN_VDestroy(sunvec_bnode)
    call FSUNMatDestroy(sunmat_Jnode)
    ierr = FSUNLinSolFree(sunls_Jnode)

    ! free items from contents, then the generic structure
    ierr = FSUNNonlinSolFree(sunnls_LOC)

    call FSUNNonlinSolFreeEmpty(sunnls)

  end function TaskLocalNewton_Free


  integer(c_int) function TaskLocalNewton_SetSysFn(sunnls, SysFn) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nonlinearsolver_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls ! nonlinear solver
    type(c_funptr)           :: SysFn  ! residual function

    !======= Internals ============

    ierr = FSUNNonlinSolSetSysFn(sunnls_LOC, SysFn)

  end function TaskLocalNewton_SetSysFn


  integer(c_int) function TaskLocalNewton_SetConvTestFn(sunnls, CTestFn, &
       ctest_data) result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nonlinearsolver_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls     ! nonlinear solver
    type(c_funptr), value    :: CTestFn    ! convergence test function
    type(c_ptr), value       :: ctest_data ! convergence test data

    !======= Internals ============

    ierr = FSUNNonlinSolSetConvTestFn(sunnls_LOC, CTestFn, ctest_data)

  end function TaskLocalNewton_SetConvTestFn


  integer(c_int) function TaskLocalNewton_GetNumConvFails(sunnls, nconvfails) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nonlinearsolver_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(SUNNonlinearSolver) :: sunnls     ! nonlinear solver
    integer(c_long)          :: nconvfails ! convergence fails

    !======= Internals ============

    nconvfails = ncnf_loc

    ! return success
    ierr = 0

  end function TaskLocalNewton_GetNumConvFails


  function TaskLocalNewton(arkode_mem, sunvec_y) result(sunnls)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod
    use fsundials_nonlinearsolver_mod
    use fsunnonlinsol_newton_mod
    use fsunlinsol_dense_mod
    use fsunmatrix_dense_mod

    use ode_mod, only : Nvar, comm, DFID, monitor

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr), target, intent(in) :: arkode_mem ! ARKODE memory
    type(N_Vector)                  :: sunvec_y   ! solution N_Vector

    type(SUNNonlinearSolver),     pointer :: sunnls ! SUNDIALS nonlinear solver
    type(SUNNonlinearSolver_Ops), pointer :: nlsops ! solver operations

    integer        :: ierr   ! MPI error status
    integer(c_int) :: retval ! SUNDIALS error status

    !======= Internals ============

    ! Set pointer to ARKODE memory structure
    arkmem => arkode_mem

    ! Create an empty nonlinear linear solver object
    sunnls => FSUNNonlinSolNewEmpty()
    if (.not. associated(sunnls)) then
       print *, "Error: FSUNNonlinSolNewEmpty returned NULL"
       call MPI_Abort(comm, 1, ierr)
    end if

    ! Access the SUNNonlinearSolver ops structure
    call c_f_pointer(sunnls%ops, nlsops)

    ! Attach operations
    nlsops%gettype         = c_funloc(TaskLocalNewton_GetType)
    nlsops%initialize      = c_funloc(TaskLocalNewton_Initialize)
    nlsops%solve           = c_funloc(TaskLocalNewton_Solve)
    nlsops%free            = c_funloc(TaskLocalNewton_Free)
    nlsops%setsysfn        = c_funloc(TaskLocalNewton_SetSysFn)
    nlsops%setctestfn      = c_funloc(TaskLocalNewton_SetConvTestFn)
    nlsops%getnumconvfails = c_funloc(TaskLocalNewton_GetNumConvFails)

    ! Create the task local Newton solver
    sunnls_LOC => FSUNNonlinSol_Newton(sunvec_y)

    ! Create vector pointers to receive residual data
    zpred_ptr = FN_VNewVectorArray(1)
    z_ptr     = FN_VNewVectorArray(1)
    Fi_ptr    = FN_VNewVectorArray(1)
    sdata_ptr = FN_VNewVectorArray(1)

    sunvec_bnode => FN_VNew_Serial(int(Nvar, c_long))
    sunmat_Jnode => FSUNDenseMatrix(int(Nvar, c_long), int(Nvar, c_long))
    sunls_Jnode  => FSUNLinSol_Dense(sunvec_bnode, sunmat_Jnode)

    ! initialize number of nonlinear solver function evals and fails
    nnlfi    = 0
    ncnf_loc = 0

    if (monitor) then
       retval = FSUNNonlinSolSetInfoFile_Newton(sunnls_LOC, DFID)
       if (retval /= 0) then
          print *, "Error: FSUNNonlinSolSetInfoFile_Newton returned ",retval
          call MPI_Abort(comm, 1, ierr)
       end if

       retval = FSUNNonlinSolSetPrintLevel_Newton(sunnls_LOC, 1)
       if (retval /= 0) then
          print *, "Error: FSUNNonlinSolSetPrintLevel_Newton returned ",retval
          call MPI_Abort(comm, 1, ierr)
       end if
    end if

  end function TaskLocalNewton

end module nls_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_types_mod        ! sundials defined types
  use fsundials_nvector_mod      ! Access generic N_Vector
  use fnvector_mpiplusx_mod      ! Access MPI+X N_Vector
  use fnvector_mpimanyvector_mod ! Access MPIManyVector N_Vector
  use fnvector_serial_mod        ! Access serial N_Vector

  use ode_mod, only : comm, myid, Nx, Neq, dx, fused, explicit, printtime, nout

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  integer          :: i
  integer          :: ierr               ! MPI error status
  integer(c_int)   :: retval             ! SUNDIALS error status
  double precision :: starttime, endtime ! timing variables

  type(N_Vector), pointer :: sunvec_ys   ! sundials serial vector
  type(N_Vector), pointer :: sunvec_y    ! sundials MPI+X vector

  !======= Internals ============

  ! Initialize MPI
  call MPI_Init(ierr)
  if (ierr /= 0) then
     print *, "Error: MPI_Init returned ",ierr
     stop 1
  end if

  ! Start timing
  starttime = MPI_Wtime()

  ! Process input args and setup the problem
  call SetupProblem()

  ! Create solution vector
  sunvec_ys => FN_VNew_Serial(int(Neq, c_long))
  sunvec_y  => FN_VMake_MPIPlusX(comm, sunvec_ys)

  ! Enable fused vector ops in local and MPI+X vectors
  if (fused) then
     retval = FN_VEnableFusedOps_Serial(sunvec_ys, SUNTRUE)
     if (retval /= 0) then
        print *, "Error: FN_VEnableFusedOps_Serial returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FN_VEnableFusedOps_MPIManyVector(sunvec_y, SUNTRUE)
     if (retval /= 0) then
        print *, "Error: FN_VEnableFusedOps_MPIManyVector returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

  ! Set the initial condition
  call SetIC(sunvec_y)

  ! Output spatial mesh to disk (add extra point for periodic BC
  if (myid == 0 .and. nout > 0) then
     open(99, file="mesh.txt")
     do i = 0, Nx
        write(99, "(es23.16)") dx * i
     end do
  end if

  ! Integrate in time
  if (explicit) then
     call EvolveProblemExplicit(sunvec_y)
  else
     call EvolveProblemIMEX(sunvec_y)
  end if

  ! End timing
  endtime = MPI_Wtime()

  if (myid == 0 .and. printtime) then
     print "(A,es12.5,A)", "Total wall clock time: ",endtime-starttime," sec"
  end if

  ! Finalize MPI
  call FN_VDestroy(sunvec_ys)
  call FN_VDestroy(sunvec_y)
  call FreeProblem()
  call MPI_Finalize(ierr)

end program main


! Setup ARKODE and evolve problem in time with IMEX method
subroutine EvolveProblemIMEX(sunvec_y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_futils_mod          ! Fortran utilities
  use farkode_mod                   ! Access ARKode
  use farkode_arkstep_mod           ! Access ARKStep
  use fsundials_nvector_mod         ! Access generic N_Vector
  use fsundials_matrix_mod          ! Access generic SUNMatrix
  use fsunmatrix_dense_mod          ! Access dense SUNMatrix
  use fsundials_linearsolver_mod    ! Access generic SUNLinearSolver
  use fsunlinsol_dense_mod          ! Access dense SUNLinearSolver
  use fsunlinsol_spgmr_mod          ! Access GMRES SUNLinearSolver
  use fsundials_nonlinearsolver_mod ! Access generic SUNNonlinearSolver
  use fsunnonlinsol_newton_mod      ! Access Newton SUNNonlinearSolver

  use ode_mod, only : comm, myid, Neq, t0, tf, atol, rtol, order, DFID, &
       monitor, global, nout, umask_s, Advection, Reaction

  use prec_mod, only : sunls_P, sunmat_P, PSetup, PSolve

  use nls_mod, only : nnlfi, TaskLocalNewton

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  ! calling variables
  type(N_Vector) :: sunvec_y ! solution N_Vector

  ! local variables
  type(c_ptr)     :: arkode_mem ! ARKODE memory
  real(c_double)  :: t(1)       ! ARKODE current time
  integer(c_int)  :: retval     ! ARKODE return value
  integer(c_long) :: nst(1)     ! number of time steps
  integer(c_long) :: nst_a(1)   ! number of step attempts
  integer(c_long) :: nfe(1)     ! number of explicit RHS evals
  integer(c_long) :: nfi(1)     ! number of implicit RHS evals
  integer(c_long) :: netf(1)    ! number of error test fails
  integer(c_long) :: nni(1)     ! number of nonlinear iters
  integer(c_long) :: ncfn(1)    ! number of convergence fails
  integer(c_long) :: nli(1)     ! number of linear iters
  integer(c_long) :: npre(1)    ! number of preconditioner setups
  integer(c_long) :: npsol(1)   ! number of preconditioner solves

  type(SUNNonlinearSolver), pointer :: sun_NLS  ! nonlinear solver
  type(SUNLinearSolver),    pointer :: sun_LS   ! linear solver
  type(SUNMatrix),          pointer :: sunmat_A ! sundials matrix

  integer            :: ierr        ! MPI error status
  integer            :: iout        ! output counter
  double precision   :: tout, dtout ! output time and update
  character(len=100) :: outname     ! diagnostics ouptput file

  !======= Internals ============

  ! Create the ARK timestepper module
  arkode_mem = FARKStepCreate(c_funloc(Advection), c_funloc(Reaction), &
       t0, sunvec_y)
  if (.not. c_associated(arkode_mem)) then
     print *, "Error: FARKStepCreate returned NULL"
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Select the method order
  retval = FARKStepSetOrder(arkode_mem, order)
  if (retval /= 0) then
     print *, "Error: FARKStepSetOrder returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify tolerances
  retval = FARKStepSStolerances(arkode_mem, rtol, atol)
  if (retval /= 0) then
     print *, "Error: FARKStepSStolerances returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Increase the max number of steps allowed between outputs
  retval = FARKStepSetMaxNumSteps(arkode_mem, int(100000, c_long))
  if (retval /= 0) then
     print *, "Error: FARKStepMaxNumSteps returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Open output file for integrator diagnostics
  if (monitor) then

     write(outname,"(A,I0.6,A)") "diagnostics.",myid,".txt"
     DFID = FSUNDIALSFileOpen(trim(outname), "w")

     retval = FARKStepSetDiagnostics(arkode_mem, DFID)
     if (retval /= 0) then
        print *, "Error: FARKStepSetDiagnostics returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

  end if

  ! Create the (non)linear solver
  if (global) then

     ! Create nonlinear solver
     sun_NLS => FSUNNonlinSol_Newton(sunvec_y)
     if (.not. associated(sun_NLS)) then
        print *, "Error: SUNNonlinSol_Newton returned NULL"
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Attach nonlinear solver
     retval = FARKStepSetNonlinearSolver(arkode_mem, sun_NLS)
     if (retval /= 0) then
        print *, "Error: FARKStepSetNonlinearSolver returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Create linear solver
     sun_LS => FSUNLinSol_SPGMR(sunvec_y, PREC_LEFT, 0)
     if (.not. associated(sun_NLS)) then
        print *, "Error: FSUNLinSol_SPGMR returned NULL"
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Attach linear solver
     sunmat_A => null()
     retval = FARKStepSetLinearSolver(arkode_mem, sun_LS, sunmat_A)
     if (retval /= 0) then
        print *, "Error: FARKStepSetLinearSolver returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Attach preconditioner
     retval = FARKStepSetPreconditioner(arkode_mem, c_funloc(PSetup), &
          c_funloc(PSolve))
     if (retval /= 0) then
        print *, "Error: FARKStepSetPreconditioner returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Create MPI task-local data structures for preconditioning
     sunmat_P => FSUNDenseMatrix(int(Neq, c_long), int(Neq, c_long))
     sunls_P  => FSUNLinSol_Dense(umask_s, sunmat_P)

  else

     ! The custom task-local nonlinear solver handles the linear solve
     ! as well, so we do not need a SUNLinearSolver
     sun_NLS => TaskLocalNewton(arkode_mem, umask_s)
     if (.not. associated(sun_NLS)) then
        print *, "Error: TaskLocalNewton returned NULL"
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Attach nonlinear solver
     retval = FARKStepSetNonlinearSolver(arkode_mem, sun_NLS)
     if (retval /= 0) then
        print *, "Error: FARKStepSetNonlinearSolver returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

  end if

  ! Set initial time, determine output time, and initialize output count
  t(1)  = t0
  dtout = (tf - t0)
  if (nout /= 0) then
     dtout = dtout / nout
  end if
  tout = t(1) + dtout
  iout = 0

  ! Output initial condition
  if (myid == 0 .and. monitor) then
     print *, ""
     print *, "     t           ||u||_rms       ||v||_rms       ||w||_rms"
     print *, "-----------------------------------------------------------"
  end if
  call WriteOutput(t, sunvec_y)

  ! Integrate to final time
  do while (iout < max(nout,1))

     ! Integrate to output time
     retval = FARKStepEvolve(arkode_mem, tout, sunvec_y, t, ARK_NORMAL)
     if (retval /= 0) then
        print *, "Error: FARKStepEvolve returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Output state
     call WriteOutput(t, sunvec_y)

     ! Update output time
     tout = tout + dtout
     if (tout > tf) then
        tout = tf
     end if

     ! Update output count
     iout = iout + 1

  end do

  if (myid == 0 .and. monitor) then
     print *, "-----------------------------------------------------------"
     print *, ""
  end if

  ! close output stream
  if (monitor) call FSUNDIALSFileClose(DFID)

  ! Get final statistics
  retval = FARKStepGetNumSteps(arkode_mem, nst)
  if (retval /= 0) then
     print *, "Error: FARKStepGetNumSteps returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKStepGetNumStepAttempts(arkode_mem, nst_a)
  if (retval /= 0) then
     print *, "Error: FARKStepGetNumStepAttempts returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKStepGetNumRhsEvals(arkode_mem, nfe, nfi)
  if (retval /= 0) then
     print *, "Error: FARKStepGetNumRhsEvals returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKStepGetNumErrTestFails(arkode_mem, netf)
  if (retval /= 0) then
     print *, "Error: FARKStepGetNumErrTestFails returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKStepGetNumNonlinSolvIters(arkode_mem, nni)
  if (retval /= 0) then
     print *, "Error: FARKStepGetNumNonlinSolvIters returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FARKStepGetNumNonlinSolvConvFails(arkode_mem, ncfn)
  if (retval /= 0) then
     print *, "Error: FARKStepGetNumNonlinSolvConvFails returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  if (global) then

     retval = FARKStepGetNumLinIters(arkode_mem, nli)
     if (retval /= 0) then
        print *, "Error: FARKStepGetNumLinIters returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FARKStepGetNumPrecEvals(arkode_mem, npre)
     if (retval /= 0) then
        print *, "Error: FARKStepGetNumPrecEvals returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FARKStepGetNumPrecSolves(arkode_mem, npsol)
     if (retval /= 0) then
        print *, "Error: FARKStepGetNumPrecSolves returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

  end if

  ! Print final statistics
  if (myid == 0) then

     print "(A)","Final Solver Statistics (for processor 0):"
     print "(2x,A,i0)", "Steps            = ",nst
     print "(2x,A,i0)", "Step attempts    = ",nst_a
     print "(2x,A,i0)", "Error test fails = ",netf
     print "(2x,A,i0)", "NLS fails        = ",ncfn

     if (global) then

        print "(2x,A,i0)", "RHS evals Fe     = ",nfe
        print "(2x,A,i0)", "RHS evals Fi     = ",nfi
        print "(2x,A,i0)", "NLS iters        = ",nni
        print "(2x,A,i0)", "LS iters         = ",nli
        print "(2x,A,i0)", "P setups         = ",npre
        print "(2x,A,i0)", "P solves         = ",npsol

     else

        print "(2x,A,i0)", "RHS evals Fe     = ",nfe
        print "(2x,A,i0)", "RHS evals Fi     = ",nfi + nnlfi

     end if

  end if

  ! Clean up
  call FARKStepFree(arkode_mem)

  ! Free nonlinear solver
  retval = FSUNNonlinSolFree(sun_NLS)
  if (retval /= 0) then
     print *, "Error: FSUNNonlinSolFree returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  if (global) then
     ! free task-local preconditioner solve structures
     call FSUNMatDestroy(sunmat_P)
     retval = FSUNLinSolFree(sunls_P)
     if (retval /= 0) then
        print *, "Error: FSUNLinSolFree returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! free global linear solver
     retval = FSUNLinSolFree(sun_LS)
     if (retval /= 0) then
        print *, "Error: FSUNLinSolFree returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

end subroutine EvolveProblemIMEX


subroutine EvolveProblemExplicit(sunvec_y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_futils_mod  ! Fortran utilities
  use farkode_mod           ! Access ARKode
  use farkode_erkstep_mod   ! Access ERKStep
  use fsundials_nvector_mod ! Access generic N_Vector

  use ode_mod, only : comm, myid, t0, tf, atol, rtol, order, DFID, monitor, &
       nout, AdvectionReaction

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  ! calling variables
  type(N_Vector) :: sunvec_y ! solution N_Vector

  ! local variables
  type(c_ptr)     :: arkode_mem ! ARKODE memory
  real(c_double)  :: t(1)       ! ARKODE current time
  integer(c_int)  :: retval     ! ARKODE return value
  integer(c_long) :: nst(1)     ! number of time step
  integer(c_long) :: nst_a(1)   ! number of step attempts
  integer(c_long) :: nfe(1)     ! number of RHS evals
  integer(c_long) :: netf(1)    ! number of error test fails

  integer            :: ierr        ! output counter
  integer            :: iout        ! output counter
  double precision   :: tout, dtout ! output time and update
  character(len=100) :: outname     ! diagnostics ouptput file

  !======= Internals ============

  ! Create the ERK integrator
  arkode_mem = FERKStepCreate(c_funloc(AdvectionReaction), t0, sunvec_y)
  if (.not. c_associated(arkode_mem)) then
     print *, "Error: FERKStepCreate returned NULL"
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Select the method order
  retval = FERKStepSetOrder(arkode_mem, order)
  if (retval /= 0) then
     print *, "Error: FERKStepSetOrder returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Specify tolerances
  retval = FERKStepSStolerances(arkode_mem, rtol, atol)
  if (retval /= 0) then
     print *, "Error: FERKStepSStolerances returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Increase the max number of steps allowed between outputs
  retval = FERKStepSetMaxNumSteps(arkode_mem, int(100000, c_long))
  if (retval /= 0) then
     print *, "Error: FERKStepMaxNumSteps returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Open output file for integrator diagnotics
  if (monitor) then

     write(outname,"(A,I0.6,A)") "diagnostics.",myid,".txt"
     DFID = FSUNDIALSFileOpen(trim(outname), "w")

     retval = FERKStepSetDiagnostics(arkode_mem, DFID)
     if (retval /= 0) then
        print *, "Error: FERKStepSetDiagnostics returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

  end if

  ! Set initial time, determine output time, and initialize output count
  t(1)  = t0
  dtout = (tf - t0)
  if (nout /= 0) then
     dtout = dtout / nout
  end if
  tout = t(1) + dtout
  iout = 0

  ! Ouput initial condition
  if (myid == 0 .and. monitor) then
     print *, ""
     print *, "     t           ||u||_rms       ||v||_rms       ||w||_rms"
     print *, "-----------------------------------------------------------"
  end if
  call WriteOutput(t, sunvec_y)

  ! Integrate to final time
  do while (iout < nout)

     ! Integrate to output time
     retval = FERKStepEvolve(arkode_mem, tout, sunvec_y, t, ARK_NORMAL)
     if (retval /= 0) then
        print *, "Error: FERKStepEvolve returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     ! Output state
     call WriteOutput(t, sunvec_y)

     ! Update output time
     tout = tout + dtout
     if (tout > tf) then
        tout = tf
     end if

     ! Update output count
     iout = iout + 1

  end do

  if (myid == 0 .and. monitor) then
     print *, "-----------------------------------------------------------"
     print *, ""
  end if

  ! close output stream
  if (monitor) call FSUNDIALSFileClose(DFID)

  ! Get final statistics
  retval = FERKStepGetNumSteps(arkode_mem, nst)
  if (retval /= 0) then
     print *, "Error: FERKStepGetNumSteps returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FERKStepGetNumStepAttempts(arkode_mem, nst_a)
  if (retval /= 0) then
     print *, "Error: FERKStepGetNumStepAttempts returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FERKStepGetNumRhsEvals(arkode_mem, nfe)
  if (retval /= 0) then
     print *, "Error: FERKStepGetNumRhsEvals returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  retval = FERKStepGetNumErrTestFails(arkode_mem, netf)
  if (retval /= 0) then
     print *, "Error: FERKStepGetNumErrTestFails returned ",retval
     call MPI_Abort(comm, 1, ierr)
  end if

  ! Print final statistcs
  if (myid == 0) then
     print "(A)","Final Solver Statistics (for processor 0):"
     print "(2x,A,i0)", "Steps            = ",nst
     print "(2x,A,i0)", "Step attempts    = ",nst_a
     print "(2x,A,i0)", "Error test fails = ",netf
     print "(2x,A,i0)", "RHS evals        = ",nfe
  end if

  ! Clean up
  call FERKStepFree(arkode_mem)

end subroutine EvolveProblemExplicit


! Write time and solution to disk
subroutine WriteOutput(t, sunvec_y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod           ! Access ARKode
  use farkode_erkstep_mod   ! Access ERKStep
  use fsundials_nvector_mod ! Access generic N_Vector

  use ode_mod, only : Nvar, nprocs, myid, Erecv, Nx, Npts, monitor,  nout, &
       umask, vmask, wmask

  !======= Declarations =========
  implicit none

  ! calling variables
  real(c_double) :: t(1)     ! current time
  type(N_Vector) :: sunvec_y ! solution N_Vector

  real(c_double), pointer :: ydata(:) ! vector data

  integer i, idx ! loop counter and array index

  double precision :: u, v, w ! RMS norm of chemical species

  !======= Internals ============

  ! output current solution norm to screen
  if (monitor) then

     u = FN_VWL2Norm(sunvec_y, umask)
     u = sqrt(u * u / Nx)

     v = FN_VWL2Norm(sunvec_y, vmask)
     v = sqrt(v * v / Nx)

     w = FN_VWL2Norm(sunvec_y, wmask)
     w = sqrt(w * w / Nx)

     if (myid == 0) then
        print "(4(es12.5,4x))", t, u, v, w
     end if

  end if

  if (nout > 0) then

     ! get left end point for output
     call ExchangeBCOnly(sunvec_y)

     ! get vector data array
     ydata => FN_VGetArrayPointer(sunvec_y)

     ! output the times to disk
     if (myid == 0) then
        write(100,"(es23.16)") t
     end if

     ! output results to disk
     do i = 1, Npts
        idx = (i - 1) * Nvar
        write(101, "(es23.16)", advance="no") ydata(idx + 1)
        write(102, "(es23.16)", advance="no") ydata(idx + 2)
        write(103, "(es23.16)", advance="no") ydata(idx + 3)
     end do

     ! we have one extra output because of the periodic BCs
     if (myid == nprocs - 1) then
        write(101,"(es23.16)") Erecv(1)
        write(102,"(es23.16)") Erecv(2)
        write(103,"(es23.16)") Erecv(3)
     else
        write(101,"(es23.16)")
        write(102,"(es23.16)")
        write(103,"(es23.16)")
     end if

  end if

end subroutine WriteOutput


! Initial Condition Function
subroutine SetIC(sunvec_y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod

  use ode_mod, only : Nvar, myid, Npts, xmax, dx, A, B, k1, k2, k4, k3

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector) :: sunvec_y ! solution N_Vector

  ! local variables
  real(c_double), pointer :: ydata(:) ! vector data

  double precision :: x, us, vs, ws       ! position and state
  double precision :: p, mu, sigma, alpha ! perturbation vars

  integer :: j, idx ! loop counter and array index

  !======= Internals ============

  ! Access data array from SUNDIALS vector
  ydata => FN_VGetArrayPointer(sunvec_y)

  ! Steady state solution
  us = k1 * A / k4
  vs = k2 * k4 * B / (k1 * k3 * A)
  ws = 3.0d0

  ! Perturbation parameters
  mu    = xmax / 2.0d0
  sigma = xmax / 4.0d0
  alpha = 0.1d0

  ! Gaussian perturbation
  do j = 1,Npts

     x = (myid * Npts + (j - 1)) * dx
     p = alpha * exp( -((x - mu) * (x - mu)) / (2.0d0 * sigma * sigma) )

     idx = (j - 1) * Nvar

     ydata(idx + 1) = us + p
     ydata(idx + 2) = vs + p
     ydata(idx + 3) = ws + p

  end do

end subroutine SetIC


! --------------------------------------------------------------
! Utility functions
! --------------------------------------------------------------


! Exchanges the periodic BCs only by sending the first mesh node to the last
! processor.
subroutine ExchangeBCOnly(sunvec_y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod

  use ode_mod, only : Nvar, comm, myid, nprocs, Erecv, Wsend

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  ! calling variables
  type(N_Vector) :: sunvec_y ! solution N_Vector

  ! local variables
  real(c_double), pointer :: ydata(:) ! array data

  integer :: stat(MPI_STATUS_SIZE) ! MPI status
  integer :: ierr, reqS, reqR      ! MPI error status and request handles
  integer :: first, last           ! MPI process IDs

  !======= Internals ============

  ! first and last MPI task ID
  first = 0
  last  = nprocs - 1

  ! Access data array from SUNDIALS vector
  ydata => FN_VGetArrayPointer(sunvec_y)

  ! open the East Irecv buffer
  if (myid == last) then
     call MPI_Irecv(Erecv, nvar, MPI_DOUBLE_PRECISION, first, MPI_ANY_TAG, &
          comm, reqR, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Irecv returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

  ! send first mesh node to the last processor
  if (myid == first) then
     Wsend(1:Nvar) = ydata(1:Nvar)
     call MPI_Isend(Wsend, nvar, MPI_DOUBLE, last, 0, &
          comm, reqS, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Isend returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

  ! wait for exchange to finish
  if (myid == last) then
     call  MPI_Wait(reqR, stat, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Wait returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

  if (myid == first) then
     call  MPI_Wait(reqS, stat, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Wait returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

end subroutine ExchangeBCOnly


! Starts the exchange of the neighbor information
subroutine ExchangeAllStart(sunvec_y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod

  use ode_mod, only : Nvar, comm, myid, nprocs, reqS, reqR, Wrecv, Wsend, &
       Erecv, Esend, Npts, c

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  ! calling variables
  type(N_Vector) :: sunvec_y ! solution N_Vector

  ! local variables
  real(c_double), pointer :: ydata(:) ! array data

  integer :: ierr        ! MPI error status
  integer :: first, last ! MPI process IDs
  integer :: ipW, ipE    ! MPI process IDs

  !======= Internals ============

  ! first and last MPI task ID
  first = 0
  last  = nprocs - 1

  ! get the ID for the process to the West and East of this process
  if (myid == first) then
     ipW = last
  else
     ipW = myid - 1
  end if

  if (myid == last) then
     ipE = first
  else
     ipE = myid + 1
  end if

  ! Access data array from SUNDIALS vector
  ydata => FN_VGetArrayPointer(sunvec_y)

  if (c > 0.0d0) then

     ! Right moving flow uses backward difference.
     ! Send from west to east (last processor sends to first)

     call MPI_Irecv(Wrecv, nvar, MPI_DOUBLE_PRECISION, ipW, &
          MPI_ANY_TAG, comm, reqR, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Irecv returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if

     Esend(1:Nvar) = ydata(Nvar * Npts - 2 : Nvar * Npts)

     call MPI_Isend(Esend, nvar, MPI_DOUBLE_PRECISION, ipE, &
          0, comm, reqS, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Isend returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if

  else if (c < 0.0d0) then

     ! Left moving flow uses forward difference.
     ! Send from east to west (first processor sends to last)

     call MPI_Irecv(Erecv, nvar, MPI_DOUBLE_PRECISION, ipE, &
          MPI_ANY_TAG, comm, reqR, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Irecv returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if

     Wsend(1:Nvar) = ydata(1:Nvar)

     call MPI_Isend(Wsend, nvar, MPI_DOUBLE_PRECISION, ipW, &
          0, comm, reqS, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Isend returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if

  end if

end subroutine ExchangeAllStart


! Completes the exchange of the neighbor information
subroutine ExchangeAllEnd()

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use ode_mod, only : comm, reqS, reqR, c

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  integer :: stat(MPI_STATUS_SIZE) ! MPI status
  integer :: ierr                  ! MPI error status

  !======= Internals ============

  ! wait for exchange to finish
  if (c < 0.0d0 .or. c > 0.0d0) then
     call MPI_Wait(reqR, stat, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Wait returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if

     call MPI_Wait(reqS, stat, ierr)
     if (ierr /= MPI_SUCCESS) then
        print *, "Error: MPI_Wait returned ",ierr
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

end subroutine ExchangeAllEnd


subroutine SetupProblem()

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_types_mod
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use fnvector_mpiplusx_mod
  use fnvector_mpimanyvector_mod

  use ode_mod

  !======= Declarations =========
  implicit none

  ! With MPI-3 use mpi_f08 is preferred
  include "mpif.h"

  ! local variables
  real(c_double), pointer :: data(:)               ! array data
  integer(c_int)          :: retval                ! SUNDIALS error status
  integer                 :: ierr                  ! MPI error status
  integer                 :: j                     ! loop counter
  integer                 :: nargs, length, status ! input parsing vars
  character(len=32)       :: arg                   ! input arg
  character(len=100)      :: outname               ! output file name

  !======= Internals ============

  ! MPI variables
  comm = MPI_COMM_WORLD

  call MPI_Comm_rank(comm, myid, ierr)
  if (ierr /= MPI_SUCCESS) then
     print *, "Error:MPI_Comm_rank = ", ierr
     call MPI_Abort(comm, 1, ierr)
  end if

  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  if (ierr /= MPI_SUCCESS) then
     print *, "Error:MPI_Comm_rank = ", ierr
     call MPI_Abort(comm, 1, ierr)
  end if

  ! default problem setting
  Nx   = 100
  Npts = Nx / nprocs
  Neq  = Nvar * Npts

  xmax = 1.0d0
  dx   = xmax / Nx

  ! Problem parameters
  c  = 0.01d0
  A  = 1.0d0
  B  = 3.5d0
  k1 = 1.0d0
  k2 = 1.0d0
  k3 = 1.0d0
  k4 = 1.0d0
  k5 = 1.0d0 / 5.0d-6
  k6 = 1.0d0 / 5.0d-6

  ! Set default integrator options
  order     = 3
  rtol      = 1.0d-6
  atol      = 1.0d-9
  t0        = 0.0d0
  tf        = 10.0d0
  explicit  = .false.
  global    = .false.
  fused     = .false.
  monitor   = .false.
  printtime = .false.
  nout      = 40

  ! check for input args
  nargs = command_argument_count()

  j = 1
  do while (j <= nargs)

     ! get input arg
     call get_command_argument(j, arg, length, status)

     ! check if reading the input was successful
     if (status == -1) then
        print *, "ERROR: Command line input too long (max length = 32)"
        call MPI_Abort(comm, 1, ierr)
     end if

     ! check if there are no more inputs to read
     if (len_trim(arg) == 0) exit

     ! check for valid input options
     if (trim(arg) == "--monitor") then
        monitor = .true.
     else if (trim(arg) == "--printtime") then
        printtime = .true.
     else if (trim(arg) == "--nout") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) nout
     else if (trim(arg) == "--nx") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) Nx
     else if (trim(arg) == "--xmax") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) xmax
     else if (trim(arg) == "--A") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) A
     else if (trim(arg) == "--B") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) B
     else if (trim(arg) == "--k") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) k1
        read(arg,*) k2
        read(arg,*) k3
        read(arg,*) k4
     else if (trim(arg) == "--c") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) c
     else if (trim(arg) == "--order") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) order
     else if (trim(arg) == "--explicit") then
        explicit = .true.
     else if (trim(arg) == "--global-nls") then
        global = .true.
     else if (trim(arg) == "--fused") then
        fused = .true.
     else if (trim(arg) == "--tf") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) tf
     else if (trim(arg) == "--rtol") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) rtol
     else if (trim(arg) == "--atol") then
        j = j + 1
        call get_command_argument(j, arg)
        read(arg,*) atol
     else if (trim(arg) == "--help") then
        if (myid == 0) call InputHelp()
        call MPI_Abort(comm, 1, ierr)
     else
        if (myid == 0) then
           print *, "Error: Unknown command line input ",trim(arg)
           call InputHelp()
        end if
        call MPI_Abort(comm, 1, ierr)
     end if

     ! move to the next input
     j = j+1
  end do

  ! Setup the parallel decomposition
  if (MOD(Nx,nprocs) > 0) then
     print *, "ERROR: The mesh size (nx = ", Nx,") must be divisible by the number of processors (",nprocs,")"
     call MPI_Abort(comm, 1, ierr)
  end if

  Npts = Nx / nprocs
  Neq  = nvar * Npts
  dx   = xmax / Nx   ! Nx is number of intervals

  ! Create the solution masks
  umask_s => FN_VNew_Serial(int(Neq, c_long))
  umask   => FN_VMake_MPIPlusX(comm, umask_s)

  if (fused) then
     retval = FN_VEnableFusedOps_Serial(umask_s, SUNTRUE)
     if (retval /= 0) then
        print *, "Error: FN_VEnableFusedOps_Serial returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if

     retval = FN_VEnableFusedOps_MPIManyVector(umask, SUNTRUE)
     if (retval /= 0) then
        print *, "Error: FN_VEnableFusedOps_MPIManyVector returned ",retval
        call MPI_Abort(comm, 1, ierr)
     end if
  end if

  call FN_VConst(0.0d0, umask)
  data => FN_VGetArrayPointer(umask)
  do j = 1, Npts
     data(1 + (j - 1) * nvar) = 1.0d0
  end do

  vmask => FN_VClone(umask)

  call FN_VConst(0.0d0, vmask)
  data => FN_VGetArrayPointer(vmask)
  do j = 1, Npts
     data(2 + (j - 1) * nvar) = 1.0d0
  end do

  wmask => FN_VClone(umask)

  call FN_VConst(0.0d0, wmask)
  data => FN_VGetArrayPointer(wmask)
  do j = 1, Npts
     data(3 + (j - 1) * nvar) = 1.0d0
  end do

  ! Open output files for results
  if (nout > 0) then

     if (myid == 0) then
        write(outname, "(A,I0.6,A)") "t.",myid,".txt"
        open(100, file=trim(outname))
     end if

     write(outname, "(A,I0.6,A)") "u.",myid,".txt"
     open(101, file=trim(outname))

     write(outname, "(A,I0.6,A)") "v.",myid,".txt"
     open(102, file=trim(outname))

     write(outname, "(A,I0.6,A)") "w.",myid,".txt"
     open(103, file=trim(outname))

  end if

  ! Print problem setup
  if (myid == 0) then

     print "(A)"       , "1D Advection-Reaction Test Problem"
     print "(A,i0)"    , "Number of Processors = ", nprocs
     print "(A)"       , "Mesh Info:"
     print "(A,i0)"    , "  Nx   = ",nx
     print "(A,i0)"    , "  Npts = ",Npts
     print "(A,es12.5)", "  xmax = ",xmax
     print "(A,es12.5)", "  dx   = ",dx
     print "(A)"       , "Problem Parameters:"
     print "(A,es12.5)", "  A = ",A
     print "(A,es12.5)", "  B = ",B
     print "(A,es12.5)", "  k = ",k1
     print "(A,es12.5)", "  c = ",c
     print "(A)"       , "Integrator Options:"
     print "(A,es12.5)", "  t0         = ", t0
     print "(A,es12.5)", "  tf         = ", tf
     print "(A,es12.5)", "  reltol     = ", rtol
     print "(A,es12.5)", "  abstol     = ", atol
     print "(A,i0)"    , "  order      = ", order
     print "(A,L1)"    , "  explicit   = ", explicit
     print "(A,L1)"    , "  fused ops  = ", fused
     if (.not. explicit) then
        print "(A,L1)","  global NLS = ", global
     end if
     print "(A,i0)"    , "  nout       = ", nout

  end if

end subroutine SetupProblem


subroutine FreeProblem()

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_nvector_mod
  use fsundials_matrix_mod
  use fsundials_linearsolver_mod

  use ode_mod, only : myid, nout, umask_s, umask, vmask, wmask

  !======= Declarations =========
  implicit none

  !======= Internals ============

  ! free solution masks
  call FN_VDestroy(umask_s)
  call FN_VDestroy(umask)
  call FN_VDestroy(vmask)
  call FN_VDestroy(wmask)

  ! close output streams
  if (nout > 0) then
     if (myid == 0) close(100)
     close(101)
     close(102)
     close(103)
  end if

end subroutine FreeProblem


subroutine InputHelp()

  print *, "Command line options:"
  print *, "  --help           prints this message"
  print *, "  --monitor        print solution information to screen (slower)"
  print *, "  --nout <int>     number of output times"
  print *, "  --explicit       use an explicit method instead of IMEX"
  print *, "  --global-nls     use a global newton nonlinear solver instead of task-local (for IMEX only)"
  print *, "  --order <int>    the method order to use"
  print *, "  --nx <int>       number of mesh points"
  print *, "  --xmax <double>  maximum value of x (size of domain)"
  print *, "  --tf <double>    final time"
  print *, "  --A <double>     A parameter value"
  print *, "  --B <double>     B parameter value"
  print *, "  --k <double>     reaction rate"
  print *, "  --c <double>     advection speed"
  print *, "  --rtol <double>  relative tolerance"
  print *, "  --atol <double>  absolute tolerance"

end subroutine InputHelp
