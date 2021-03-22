! -----------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
! -----------------------------------------------------------------
! Based on kinLaplace_picard_bnd.c by Carol S. Woodward @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------
! This example solves a 2D elliptic PDE
!
!    d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u - 2.0
!
! subject to homogeneous Dirichlet boundary conditions.
! The PDE is discretized on a uniform NX+2 by NY+2 grid with
! central differencing, and with boundary values eliminated,
! leaving a system of size NEQ = NX*NY.
! The nonlinear system is solved by KINSOL using the Picard
! iteration and the SPGMR linear solver.
! -----------------------------------------------------------------

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
  ! Nonlinear residual function
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


  ! ----------------------------------------------------------------
  ! Jacobian vector product function
  ! ----------------------------------------------------------------

  integer(c_int) function jactimes(sunvec_v, sunvec_Jv, sunvec_u, new_u, user_data) &
       result(ierr) bind(C)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    type(N_Vector)       :: sunvec_v   ! input vector
    type(N_Vector)       :: sunvec_Jv  ! output vector
    type(N_Vector)       :: sunvec_u   ! current solution vector
    integer(c_int)       :: new_u      ! flag indicating if u has been updated
    type(c_ptr),   value :: user_data  ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: v(:,:), Jv(:,:)

    ! internal variables
    integer(c_long) :: i, j
    real(c_double)  :: dx, dy, hdiff, vdiff, hdc, vdc, vij, vdn, vup, vlt, vrt

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors, casting as 2D Fortran arrays
    v(1:nx, 1:ny)  => FN_VGetArrayPointer(sunvec_v)
    Jv(1:nx, 1:ny) => FN_VGetArrayPointer(sunvec_Jv)

    ! set shortcut constants
    dx = 1.d0/(nx+1)
    dy = 1.d0/(ny+1)
    hdc = 1.d0/(dx*dx)
    vdc = 1.d0/(dy*dy)

    ! loop over domain, computing residual
    do j = 1,ny
       do i = 1,nx

          ! Extract v at x_i, y_j and four neighboring points
          vij = v(i,j)
          vdn = 0.d0
          if (j > 1)  vdn = v(i,j-1)
          vup = 0.d0
          if (j < ny) vup = v(i,j+1)
          vlt = 0.d0
          if (i > 1)  vlt = v(i-1,j)
          vrt = 0.d0
          if (i < nx) vrt = v(i+1,j)

          ! Evaluate diffusion components
          hdiff = hdc*(vlt - 2.d0*vij + vrt)
          vdiff = vdc*(vup - 2.d0*vij + vdn)

          ! Set residual at x_i, y_j
          Jv(i, j) = hdiff + vdiff

       end do
    end do

    ! return success
    ierr = 0
    return

  end function jactimes

end module prob_mod


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fsundials_futils_mod       ! Fortran utilities
  use fkinsol_mod                ! Fortran interface to KINSOL
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunlinsol_spgmr_mod       ! Fortran interface to SPGMR SUNLinearSolver
  use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
  use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use prob_mod                   ! problem-defining functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double)  :: fnormtol, fnorm(1)
  integer(c_int)  :: ierr
  integer(c_long) :: maa = 3

  type(N_Vector),        pointer :: sunvec_u      ! sundials vectors
  type(N_Vector),        pointer :: sunvec_s
  type(SUNMatrix),       pointer :: sunmat_L      ! sundials matrix (empty)
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver

  type(c_ptr) :: kmem   ! KINSOL memory
  type(c_ptr) :: infofp ! info file

  ! solution and scaling vectors; nx, ny are set in the prob_mod module
  real(c_double), dimension(nx,ny) :: u, scale

  !======= Internals ============

  ! -------------------------
  ! Print problem description

  print *, " "
  print '(A)', "2D elliptic PDE on unit square"
  print '(A)', "   d^2 u / dx^2 + d^2 u / dy^2 = u^3 - u + 2.0"
  print '(A)', " + homogeneous Dirichlet boundary conditions"
  print *, " "
  print '(A)', "Solution method: Anderson accelerated Picard iteration with SPGMR linear solver."
  print '(2(a,i2),a,i4)', "Problem size: ", nx, " x ", ny, " = ", neq

  ! -------------------------
  ! Set initial guess, and disable scaling

  u = 0.d0
  u(2,2) = 1.0d0

  scale = 1.d0 ! no scaling used

  ! -------------------------
  ! Create vectors for solution and scaling

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

  ! Use acceleration with up to 3 prior residuals
  ierr = FKINSetMAA(kmem, maa);
  if (ierr /= 0) then
     print *, 'Error in FKINISetMAA, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINInit(kmem, c_funloc(func), sunvec_u)
  if (ierr /= 0) then
     print *, 'Error in FKINInit, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set optional inputs

  ! Specify stopping tolerance based on residual

  fnormtol = ftol
  ierr = FKINSetFuncNormTol(kmem, fnormtol)
  if (ierr /= 0) then
     print *, 'Error in FKINSetFuncNormTol, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! Set information file

  infofp = FSUNDIALSFileOpen("KINSOL.log", "w");

  ierr = FKINSetInfoFile(kmem, infofp);
  if (ierr /= 0) then
     print *, 'Error in FKINSetInfoFile, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINSetPrintLevel(kmem, 3);
  if (ierr /= 0) then
     print *, 'Error in FKINSetPrintLevel, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Create a linear solver

  sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_u, PREC_NONE, 10)
  if (.not. associated(sunlinsol_LS)) then
     print *,'ERROR: sunlinsol = NULL'
     stop 1
  end if

  ! -------------------------
  ! Attach linear solver

  sunmat_L => null()

  ierr = FKINSetLinearSolver(kmem, sunlinsol_LS, sunmat_L)
  if (ierr /= 0) then
     print *, 'Error in FKINSetLinearSolver, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! -------------------------
  ! Set Jacobian vector product function

  ierr = FKINSetJacTimesVecFn(kmem, c_funloc(jactimes));
  if (ierr /= 0) then
     print *, 'Error in FKINSetJacTimesVecFn, ierr = ', ierr, '; halting'
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

  ierr = FKINSol(kmem, sunvec_u, KIN_PICARD, sunvec_s, sunvec_s)
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
  print '(A,ES11.5,A)', "Computed solution (||F|| = ", fnorm,"):"
  print *, " "
  call PrintOutput(u)
  call PrintFinalStats(kmem)

  ! clean up
  call FSUNDIALSFileClose(infofp)
  call FKINFree(kmem)
  ierr = FSUNLinSolFree(sunlinsol_LS)
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

  write(*,'(11x)',advance='no')
  do i = 1,nx,skip
     x = i*dx
     write(*,'(f8.5,1x)',advance='no') x
  end do
  print *, " "
  print *, " "

  do j = 1,ny,skip
     y = j*dy
     write(*,'(f7.5,4x)',advance='no') y
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
  integer(c_long) :: nni(1), nli(1), ncfl(1), nfe(1), nfeLS(1), njvevals(1)
  integer(c_long) :: npe(1), nps(1), lenrw(1), leniw(1), lenrwLS(1), leniwLS(1)

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

  ! Linear solver statistics

  ierr = FKINGetNumLinIters(kmem, nli)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumLinFuncEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumLinFuncEvals(kmem, nfeLS)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumLinFuncEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumLinConvFails(kmem, ncfl)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumLinConvFails, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumJtimesEvals(kmem, njvevals)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumJtimesEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumPrecEvals(kmem, npe)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumPrecEvals, ierr = ', ierr, '; halting'
     stop 1
  end if

  ierr = FKINGetNumPrecSolves(kmem, nps)
  if (ierr /= 0) then
     print *, 'Error in FKINGetNumPrecSolves, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! Main solver workspace size

  ierr = FKINGetWorkSpace(kmem, lenrw, leniw)
  if (ierr /= 0) then
     print *, 'Error in FKINGetWorkSpace, ierr = ', ierr, '; halting'
     stop 1
  end if

  ! Linear solver workspace size

  ierr = FKINGetLinWorkSpace(kmem, lenrwLS, leniwLS)
  if (ierr /= 0) then
     print *, 'Error in FKINGetLinWorkSpace, ierr = ', ierr, '; halting'
     stop 1
  end if

  print *, ' '
  print '(A)', 'Final Statistics..'
  print *, ' '
  print '(3(A,i6))'    ,'nni = ', nni,      '  nli   = ', nli,      '  ncfl = ',ncfl
  print '(3(A,i6))'    ,'nfe = ', nfe,      '  nfeLS = ', nfeLS,    '  njt  = ',njvevals
  print '(2(A,i6))'    ,'npe = ', npe,      '  nps   = ', nps
  print *, ' '
  print '(2(A,i6))'    ,'lenrw   = ', lenrw,    '  leniw   = ', leniw
  print '(2(A,i6))'    ,'lenrwLS = ', lenrwLS,  '  leniwLS = ', leniwLS

  return

end subroutine PrintFinalStats
