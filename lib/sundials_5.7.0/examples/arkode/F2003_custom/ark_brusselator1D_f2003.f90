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
! Example problem:
!
! The following test simulates a brusselator problem from chemical
! kinetics.  This is a PDE system with 3 components, Y = [u,v,w],
! satisfying the equations,
!    u_t = du*u_xx + a - (w+1)*u + v*u^2
!    v_t = dv*v_xx + w*u - v*u^2
!    w_t = dw*w_xx + (b-w)/ep - w*u
! for t in [0, 80], x in [0, 1], with initial conditions
!    u(0,x) =  a  + 0.1*sin(pi*x)
!    v(0,x) = b/a + 0.1*sin(pi*x)
!    w(0,x) =  b  + 0.1*sin(pi*x),
! and with stationary boundary conditions, i.e.
!    u_t(t,0) = u_t(t,1) = 0,
!    v_t(t,0) = v_t(t,1) = 0,
!    w_t(t,0) = w_t(t,1) = 0.
!
! The spatial derivatives are computed using second-order
! centered differences, with the data distributed over N points
! on a uniform spatial grid.
!
! This program solves the problem with the DIRK method, using a
! Newton iteration with a custom linear solver module, and that
! uses the custom fnvector_fortran_mod vector module.
!
! 100 outputs are printed at equal intervals, and run statistics
! are printed at the end.
! ------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none
  integer(c_long), parameter :: N = 201         ! number of intervals
  integer(c_long), parameter :: Nt = 100        ! total number of output times
  integer(c_long), parameter :: Nvar = 3        ! number of solution fields
  integer(c_long), parameter :: neq = N*Nvar    ! total size of solution vector
  real(c_double),  parameter :: dx = 1.d0/(N-1) ! mesh spacing
  real(c_double),  parameter :: a = 0.6d0       ! constant forcing on u
  real(c_double),  parameter :: b = 2.d0        ! steady-state value of w
  real(c_double),  parameter :: du = 0.001d0    ! diffusion coeff for u
  real(c_double),  parameter :: dv = 0.001d0    ! diffusion coeff for v
  real(c_double),  parameter :: dw = 0.001d0    ! diffusion coeff for w
  real(c_double),  parameter :: ep = 1.d-5      ! stiffness parameter
  real(c_double),  parameter :: T0 = 0.d0       ! initial time
  real(c_double),  parameter :: Tf = 10.d0      ! final time
  real(c_double),  parameter :: reltol = 1.d-6  ! solver tolerances
  real(c_double),  parameter :: abstol = 1.d-10

contains

  ! ----------------------------------------------------------------
  ! RhsImplicit provides the right hand side function for the
  ! implicit portion of the ODE: dy/dt = fi(t,y) + fe(t,y).  Here,
  ! this contains only the 'reaction' components of the problem.
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function RhsImplicit(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='RhsImplicit')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_fortran_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn  ! current time
    type(N_Vector) :: sunvec_y   ! solution N_Vector
    type(N_Vector) :: sunvec_f   ! rhs N_Vector
    type(c_ptr)    :: user_data  ! user-defined data

    ! local variables
    type(FVec), pointer :: y, f  ! ptrs to Fortran vector data
    real(c_double)      :: u, v, w
    integer(c_long)     :: i

    !======= Internals ============

    ! extract Fortran vector structures to work with
    call c_f_pointer(sunvec_y%content, y)
    call c_f_pointer(sunvec_f%content, f)

    ! iterate over domain, computing all equations
    do i = 2,N-1

       ! set solution variable shortcuts
       u  = y%data(1,i)
       v  = y%data(2,i)
       w  = y%data(3,i)

       ! Fill in ODE RHS for u
       f%data(1,i) = a - (w+1.d0)*u + v*u*u

       ! Fill in ODE RHS for v
       f%data(2,i) = w*u - v*u*u

       ! Fill in ODE RHS for w
       f%data(3,i) = (b-w)/ep - w*u

    end do

    ! enforce stationary boundary conditions
    f%data(:,1) = 0.d0
    f%data(:,N) = 0.d0

    ! return success
    ierr = 0

  end function RhsImplicit

  ! ----------------------------------------------------------------
  ! RhsExplicit provides the right hand side function for the
  ! explicit portion of the ODE: dy/dt = fi(t,y) + fe(t,y).  Here,
  ! this contains only the 'diffusion' components of the problem.
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function RhsExplicit(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='RhsExplicit')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_fortran_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn  ! current time
    type(N_Vector) :: sunvec_y   ! solution N_Vector
    type(N_Vector) :: sunvec_f   ! rhs N_Vector
    type(c_ptr) :: user_data     ! user-defined data

    ! local variables
    type(FVec), pointer :: y, f  ! ptrs to Fortran vector data
    real(c_double)  :: dconst(3)
    integer(c_long) :: i, j

    !======= Internals ============

    ! extract Fortran vector structures to work with
    call c_f_pointer(sunvec_y%content, y)
    call c_f_pointer(sunvec_f%content, f)

    ! set shortcut variables
    dconst = (/ du/dx/dx, dv/dx/dx, dw/dx/dx /)

    ! Fill in diffusion RHS for each species
    do j = 2,N-1
       do i = 1,Nvar
          f%data(i,j) = (y%data(i,j-1) + y%data(i,j+1) - 2.d0*y%data(i,j))*dconst(i)
       end do
    end do

    ! enforce stationary boundary conditions
    f%data(:,1) = 0.d0
    f%data(:,N) = 0.d0

    ! return success
    ierr = 0

  end function RhsExplicit

  ! ----------------------------------------------------------------
  ! JacFn provides the Jacobian construction function for the
  ! implicit portion of the ODE: J = dfi/dy
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function JacFn(tn, sunvec_y, sunvec_f, sunmat_J, &
       user_data, sunvec_t1, sunvec_t2, sunvec_t3) &
       result(ierr) bind(C,name='JacFn')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsunmatrix_fortran_mod
    use fnvector_fortran_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn  ! current time
    type(N_Vector)  :: sunvec_y  ! solution N_Vector
    type(N_Vector)  :: sunvec_f  ! rhs N_Vector
    type(SUNMatrix) :: sunmat_J  ! Jacobian matrix
    type(c_ptr)     :: user_data ! user-defined data
    type(N_Vector)  :: sunvec_t1  ! temporary N_Vectors
    type(N_Vector)  :: sunvec_t2
    type(N_Vector)  :: sunvec_t3

    ! local variables
    type(FVec), pointer :: y, f  ! ptrs to Fortran vector data
    type(FMat), pointer :: J     ! ptr to Fortran matrix data
    real(c_double)      :: u, v, w
    integer(c_long)     :: i

    !======= Internals ============

    ! extract Fortran vector and matrix structures to work with
    call c_f_pointer(sunvec_y%content, y)
    call c_f_pointer(sunvec_f%content, f)
    call c_f_pointer(sunmat_J%content, J)

    ! iterate over domain, computing all Jacobian entries
    do i = 2,N-1

       ! set solution variable shortcuts
       u  = y%data(1,i)
       v  = y%data(2,i)
       w  = y%data(3,i)

       ! Fill in Jacobian of all variables wrt u
       J%data(:,1,i) = (/ 2.d0*v*u - (w+1.d0), w - 2.d0*v*u, -w /)

       ! Fill in Jacobian of all variables wrt v
       J%data(:,2,i) = (/ u*u, -u*u, 0.d0 /)

       ! Fill in Jacobian of all variables wrt w
       J%data(:,3,i) = (/ -u, u, -1.d0/ep - u/)

    end do

    ! stationary boundary conditions
    J%data(:,:,1) = 0.d0
    J%data(:,:,N) = 0.d0

    ! return success
    ierr = 0

  end function JacFn

end module ode_mod

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod                ! Fortran interface to the ARKode module
  use farkode_arkstep_mod        ! Fortran interface to the ARKStep module
  use fnvector_fortran_mod       ! Custom Fortran N_Vector
  use fsunmatrix_fortran_mod     ! Custom Fortran SUNMatrix
  use fsunlinsol_fortran_mod     ! Custom Fortran SUNLinearSolver
  use ode_mod                    ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  integer(c_int) :: ierr, iout, i
  real(c_double) :: pi, tcur(1), dTout, tout

  type(N_Vector),        pointer :: sunvec_y      ! sundials vector
  type(SUNMatrix),       pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunls         ! sundials linear solver
  type(c_ptr)                    :: arkode_mem    ! ARKODE memory

  ! solution vector, N and Nvar are set in the ode_mod moduel
  real(c_double) :: y(Nvar,N)

  !======= Internals ============

  ! initial problem output
  print *, "  "
  print *, "1D Brusselator PDE test problem (F2003 with custom data structures):"
  print '(a,i5,a,i5)', "    N = ", N, ",  NEQ = ", neq
  print '(3(a,es9.2))', "    problem parameters:  a = ", a, ",  b = ", b, ",  ep = ", ep
  print '(3(a,es10.3))', "    diffusion coefficients:  du = ",du,",  dv = ",dv,",  dw = ",dw
  print '(2(a,es8.1))', "    reltol = ",reltol,",  abstol = ",abstol

  ! initialize SUNDIALS solution vector
  sunvec_y => FN_VMake_Fortran(Nvar, N, y)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  ! Set initial conditions into y
  pi = 4.d0*datan(1.d0)
  do i = 1,N
     y(1,i) =  a  + 0.1d0*sin(pi*i*dx)  ! u
     y(2,i) = b/a + 0.1d0*sin(pi*i*dx)  ! v
     y(3,i) =  b  + 0.1d0*sin(pi*i*dx)  ! w
  end do

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_funloc(RhsExplicit), c_funloc(RhsImplicit), T0, sunvec_y)
  if (.not. c_associated(arkode_mem)) then
     print *,'ERROR: arkode_mem = NULL'
     stop 1
  end if

  ! set routines
  ierr = FARKStepSStolerances(arkode_mem, reltol, abstol)
  if (ierr /= 0) then
    write(*,*) 'Error in FARKStepSStolerances'
    stop 1
  end if

  ! initialize custom matrix data structure and solver; attach to ARKStep
  sunmat_A => FSUNMatNew_Fortran(Nvar, N)
  if (.not. associated(sunmat_A)) then
     print *,'ERROR: sunmat = NULL'
     stop 1
  end if

  sunls => FSUNLinSolNew_Fortran(Nvar, N)
  if (.not. associated(sunls)) then
     print *,'ERROR: sunls = NULL'
     stop 1
  end if

  ! Attach matrix, linear solver, and Jacobian routine to linear solver interface
  ierr = FARKStepSetLinearSolver(arkode_mem, sunls, sunmat_A)
  if (ierr /= 0) then
    write(*,*) 'Error in FARKStepSetLinearSolver'
    stop 1
  end if

  ierr = FARKStepSetJacFn(arkode_mem, c_funloc(JacFn))
  if (ierr /= 0) then
    write(*,*) 'Error in FARKStepSetJacFn'
    stop 1
  end if

  ! main time-stepping loop: calls FARKStepEvolve to perform the integration, then
  ! prints results.  Stops when the final time has been reached
  tcur(1) = T0
  dTout = (Tf-T0)/Nt
  tout = T0+dTout
  print *, "        t      ||u||_rms   ||v||_rms   ||w||_rms"
  print *, "   ----------------------------------------------"
  do iout = 1,Nt

     ! call integrator
     ierr = FARKStepEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
     if (ierr /= 0) then
        write(*,*) 'Error in FARKStepEvolve, ierr = ', ierr, '; halting'
        stop 1
     endif

     ! access/print solution statistics
     print '(4(2x,f10.6))', tcur(1), dsqrt(sum(y(1,:)**2)/N), &
          dsqrt(sum(y(2,:)**2)/N), dsqrt(sum(y(3,:)**2)/N)

     ! update output time
     tout = min(tout + dTout, Tf)

  end do
  print *, "   ----------------------------------------------"

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
! Print ARKODE statstics to standard out
! ----------------------------------------------------------------
subroutine ARKStepStats(arkode_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use farkode_mod
  use farkode_arkstep_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: arkode_mem ! solver memory structure

  integer(c_int)  :: ierr ! error flag

  integer(c_long) :: nsteps(1)     ! num steps
  integer(c_long) :: nst_a(1)      ! num steps attempted
  integer(c_long) :: nfe(1)        ! num explicit function evals
  integer(c_long) :: nfi(1)        ! num implicit function evals
  integer(c_long) :: nfeLS(1)      ! num RHS evals to setup linear solver
  integer(c_long) :: nlinsetups(1) ! num linear solver setups
  integer(c_long) :: nje(1)        ! num Jacobian evaluations
  integer(c_long) :: netfails(1)   ! num error test fails
  integer(c_long) :: nniters(1)    ! nonlinear solver iterations
  integer(c_long) :: nncfails(1)   ! nonlinear solver fails
  integer(c_long) :: njacevals(1)  ! number of Jacobian evaluations

  !======= Internals ============

  ierr = FARKStepGetNumSteps(arkode_mem, nsteps)
  ierr = FARKStepGetNumStepAttempts(arkode_mem, nst_a)
  ierr = FARKStepGetNumRhsEvals(arkode_mem, nfe, nfi)
  ierr = FARKStepGetNumLinRhsEvals(arkode_mem, nfeLS)
  ierr = FARKStepGetNumLinSolvSetups(arkode_mem, nlinsetups)
  ierr = FARKStepGetNumJacEvals(arkode_mem, nje)
  ierr = FARKStepGetNumErrTestFails(arkode_mem, netfails)
  ierr = FARKStepGetNumNonlinSolvIters(arkode_mem, nniters)
  ierr = FARKStepGetNumNonlinSolvConvFails(arkode_mem, nncfails)
  ierr = FARKStepGetNumJacEvals(arkode_mem, njacevals)

  print *, ' '
  print *, 'Final Solver Statistics:'
  print '(4x,2(A,i9),A)' ,'Internal solver steps = ',nsteps(1),', (attempted = ',nst_a(1),')'
  print '(4x,2(A,i9))'   ,'Total RHS evals: Fe = ',nfe(1),',  Fi = ',nfi(1)
  print '(4x,A,i9)'      ,'Total linear solver setups    =',nlinsetups(1)
  print '(4x,A,i9)'      ,'Total RHS evals for setting up the linear system =',nfeLS(1)
  print '(4x,A,i9)'      ,'Total number of Jacobian evaluations =',nje(1)
  print '(4x,A,i9)'      ,'Total number of Newton iterations =',nniters(1)
  print '(4x,A,i9)'      ,'Total number of nonlinear solver convergence failures =',nncfails(1)
  print '(4x,A,i9)'      ,'Total number of error test failures =',netfails(1)
  print *, ' '

  return

end subroutine ARKStepStats
