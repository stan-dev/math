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
! This example solves a 2D elliptic PDE
!
!    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
!
! subject to homogeneous Dirichlet boundary conditions.
! The PDE is discretized on a uniform NX+2 by NY+2 grid with
! central differencing, and with boundary values eliminated,
! leaving a system of size NEQ = NX*NY.
! The nonlinear system is solved by KINSOL using the SUNBAND linear
! solver.
! ------------------------------------------------------------------

module prob_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  integer(c_long), parameter :: nx = 31
  integer(c_long), parameter :: ny = 31
  integer(c_long), parameter :: neq = nx*ny
  integer(c_long), parameter :: skip = 3
  real(c_double),  parameter :: ftol = 1.d-12

contains

  ! ----------------------------------------------------------------
  ! func: The nonlinear residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function func(sunvec_u, sunvec_f, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_u  ! solution N_Vector
    type(N_Vector)       :: sunvec_f  ! rhs N_Vector
    type(c_ptr),   value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: u(:,:), f(:,:)

    ! internal variables
    integer(c_long) :: i, j
    real(c_double)  :: dx, dy, hdiff, vdiff, hdc, vdc, uij, udn, uup, ult, urt

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors, casting as 2D Fortran arrays
    u(1:nx, 1:ny) => FN_VGetArrayPointer(sunvec_u)
    f(1:nx, 1:ny) => FN_VGetArrayPointer(sunvec_f)

    ! set shortcut constants
    dx = 1.d0/(nx+1)
    dy = 1.d0/(ny+1)
    hdc = 1.d0/(dx*dx)
    vdc = 1.d0/(dy*dy)

    ! loop over domain, computing residual
    do j = 1,ny
       do i = 1,nx

          ! Extract u at x_i, y_j and four neighboring points
          uij = u(i,j)
          udn = 0.d0
          if (j > 1)  udn = u(i,j-1)
          uup = 0.d0
          if (j < ny) uup = u(i,j+1)
          ult = 0.d0
          if (i > 1)  ult = u(i-1,j)
          urt = 0.d0
          if (i < nx) urt = u(i+1,j)

          ! Evaluate diffusion components
          hdiff = hdc*(ult - 2.d0*uij + urt)
          vdiff = vdc*(uup - 2.d0*uij + udn)

          ! Set residual at x_i, y_j
          f(i, j) = hdiff + vdiff + uij - uij*uij*uij + 2.d0
       end do
    end do

    ! return success
    ierr = 0
    return

  end function func

end module prob_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fkinsol_mod                ! Fortran interface to KINSOL
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_band_mod        ! Fortran interface to band SUNMatrix
  use fsunlinsol_band_mod        ! Fortran interface to band SUNLinearSolver
  use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use prob_mod                   ! problem-defining functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double)  :: fnormtol, fnorm(1)
  integer(c_int)  :: ierr
  integer(c_long) :: mset, msubset

  type(N_Vector),        pointer :: sunvec_u      ! sundials vectors
  type(N_Vector),        pointer :: sunvec_s
  type(SUNMatrix),       pointer :: sunmat_J      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver

  type(c_ptr) :: kmem ! KINSOL memory

  ! solution and scaling vectors; nx, ny are set in the prob_mod module
  real(c_double), dimension(nx,ny) :: u, scale

  !======= Internals ============

  ! -------------------------
  ! Print problem description

  print *, " "
  print *, "2D elliptic PDE on unit square"
  print *, "   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0"
  print *, " + homogeneous Dirichlet boundary conditions"
  print *, " "
  print *, "Solution method: Modified Newton with band linear solver"
  print '(2(a,i2),a,i3)', "Problem size: ", nx, " x ", ny, " = ", neq

  ! -------------------------
  ! Set initial guess, and disable scaling

  u = 0.d0
  scale = 1.d0

  ! -------------------------
  ! Create vectors for solution and scales

  sunvec_u => FN_VMake_Serial(neq, u)
  if (.not. associated(sunvec_u)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_s => FN_VMake_Serial(neq, scale)
  if (.not. associated(sunvec_s)) then
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

  ! sunvec_u is used as a template

  ierr = FKINInit(kmem, c_funloc(func), sunvec_u)
  if (ierr /= 0) then
     print *, 'Error in FKINInit, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set optional inputs

  fnormtol = ftol
  ierr = FKINSetFuncNormTol(kmem, fnormtol)
  if (ierr /= 0) then
     print *, 'Error in FKINSetFuncNormTol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Create a band matrix

  sunmat_J => FSUNBandMatrix(neq, nx, nx)
  if (.not. associated(sunmat_J)) then
     print *,'ERROR: sunmat = NULL'
     stop 1
  end if

  ! -------------------------
  ! Create a band linear solver

  sunlinsol_LS => FSUNLinSol_Band(sunvec_u, sunmat_J)
  if (.not. associated(sunlinsol_LS)) then
     print *,'ERROR: sunlinsol = NULL'
     stop 1
  end if

  ! -------------------------
  ! Attach band linear solver

  ierr = FKINSetLinearSolver(kmem, sunlinsol_LS, sunmat_J)
  if (ierr /= 0) then
     print *, 'Error in FKINSetLinearSolver, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Parameters for Modified Newton

  ! Force a Jacobian re-evaluation every mset iterations
  mset = 100
  ierr = FKINSetMaxSetupCalls(kmem, mset)
  if (ierr /= 0) then
     print *, 'Error in FKINSetMaxSetupCalls, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! Every msubset iterations, test if a Jacobian evaluation is necessary
  msubset = 1
  ierr = FKINSetMaxSubSetupCalls(kmem, msubset)
  if (ierr /= 0) then
     print *, 'Error in FKINSetMaxSubSetupCalls, ierr = ', ierr, '; halting'
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

  ierr = FKINSol(kmem, sunvec_u, KIN_LINESEARCH, sunvec_s, sunvec_s)
  if (ierr /= 0) then
     print *, 'Error in FKINSol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Print solution and solver statistics

  ! Get scaled norm of the system function
  ierr = FKINGetFuncNorm(kmem, fnorm)
  if (ierr /= 0) then
     print *, 'Error in FKINGetFuncNorm, ierr = ', ierr, '; halting'
     stop 1
  end if
  print *, " "
  print *, "Computed solution (||F|| = ", fnorm,"):"
  print *, " "
  call PrintOutput(u)
  call PrintFinalStats(kmem)

  ! clean up
  call FKINFree(kmem)
  ierr = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_J)
  call FN_VDestroy(sunvec_u)
  call FN_VDestroy(sunvec_s)

end program main


! ----------------------------------------------------------------
! PrintOutput: prints solution at selected points
! ----------------------------------------------------------------
subroutine PrintOutput(u)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use prob_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double), dimension(nx,ny) :: u

  ! internal variables
  integer(c_long) :: i, j
  real(c_double)  :: dx, dy, x, y

  !======= Internals ============

  ! set shortcuts
  dx = 1.d0/(nx+1)
  dy = 1.d0/(ny+1)

  write(*,'(13x)',advance='no')
  do i = 1,nx,skip
     x = i*dx
     write(*,'(f8.5,1x)',advance='no') x
  end do
  print *, " "
  print *, " "

  do j = 1,ny,skip
     y = j*dy
     write(*,'(f8.5,5x)',advance='no') y
     do i = 1,nx,skip
        write(*,'(f8.5,1x)',advance='no') u(i,j)
     end do
     print *, " "
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
  integer(c_long) :: nni(1), nfe(1), nje(1), nfeB(1), lenrw(1), leniw(1)
  integer(c_long) :: lenrwB(1), leniwB(1), nbcfails(1), nbacktr(1)

  !======= Internals ============

  ! Main solver statistics

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


  ! Linesearch statistics

  ierr = FKINGetNumBetaCondFails(kmem, nbcfails)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumBetaCondFails, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumBacktrackOps(kmem, nbacktr)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumBacktrackOps, ierr = ', ierr, '; halting'
     stop 1
  end if


  ! Main solver workspace size

  ierr = FKINGetWorkSpace(kmem, lenrw, leniw)
  if (ierr /= 0) then
     print *, 'Error in FKINGetWorkSpace, ierr = ', ierr, '; halting'
     stop 1
  end if


  ! Band linear solver statistics

  ierr = FKINGetNumJacEvals(kmem, nje)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumJacEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumLinFuncEvals(kmem, nfeB)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumLinFuncEvals, ierr = ', ierr, '; halting'
     stop 1
  end if


  ! Band linear solver workspace size

  ierr = FKINGetLinWorkSpace(kmem, lenrwB, leniwB)
  if (ierr /= 0) then
     print *, 'Error in FKINGetLinWorkSpace, ierr = ', ierr, '; halting'
     stop 1
  end if

  print *, ' '
  print *, 'Final Statistics..'
  print *, ' '
  print '(2(A,i6))'    ,'nni      =', nni,      '    nfe     =', nfe
  print '(2(A,i6))'    ,'nbcfails =', nbcfails, '    nbacktr =', nbacktr
  print '(2(A,i6))'    ,'nje      =', nje,      '    nfeB    =', nfeB
  print *, ' '
  print '(2(A,i6))'    ,'lenrw    =', lenrw,    '    leniw   =', leniw
  print '(2(A,i6))'    ,'lenrwB   =', lenrwB,   '    leniwB  =', leniwB

  return

end subroutine PrintFinalStats
