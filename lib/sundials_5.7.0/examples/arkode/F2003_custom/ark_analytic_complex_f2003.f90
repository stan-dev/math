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
! The following is a simple example problem with analytical
! solution,
!    dy/dt = lamda*y + 1/(1+t^2) - lamda*atan(t)
!
! The stiffness of the problem is directly proportional to the
! value of "lamda", which is specified through an input file.  The
! real part of lamda should be negative to result in a well-posed
! ODE; for lambdas with magnitude larger than 100 the problem
! becomes quite stiff.
!
! Here we choose lamda = -0.1 + 10 i.
!
! This program solves the problem with the ERK method and
! a user-supplied complex vector module.
! Output is printed every 1.0 units of time (10 total).
! Run statistics (optional outputs) are printed at the end.
! ------------------------------------------------------------------

module ode_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none
  integer(c_long),            parameter :: neq = 1
  integer(c_int),             parameter :: Nt = 10
  complex(c_double_complex),  parameter :: lambda = (-1d-2, 10.d0)
  real(c_double),             parameter :: T0 = 0.d0
  real(c_double),             parameter :: Tf = 10.d0
  real(c_double),             parameter :: dtmax = 0.01d0
  real(c_double),             parameter :: reltol = 1.d-6
  real(c_double),             parameter :: abstol = 1.d-10

contains

  ! ----------------------------------------------------------------
  ! Rhs provides the right hand side function for the ODE:
  ! dy/dt = f(t,y).
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function Rhs(tn, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='Rhs')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_complex_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn  ! current time
    type(N_Vector) :: sunvec_y   ! solution N_Vector
    type(N_Vector) :: sunvec_f   ! rhs N_Vector
    type(c_ptr)    :: user_data  ! user-defined data

    ! local variables
    type(FVec), pointer :: y, f  ! ptrs to Fortran vector data

    !======= Internals ============

    ! extract Fortran vector structures to work with
    y => FN_VGetFVec(sunvec_y)
    f => FN_VGetFVec(sunvec_f)

    ! evaluate ODE RHS and return success
    f%data(1) = lambda*y%data(1)
    ierr = 0

  end function Rhs

  ! ----------------------------------------------------------------
  ! Sol provides the analytical solution to the ODE.
  ! ----------------------------------------------------------------
  complex(c_double_complex) function Sol(tn) result(y)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fnvector_complex_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tn  ! current time

    !======= Internals ============

    ! fill analytical solution for this value of tn, and return
    y = 2.d0*exp(lambda*tn)
    return

  end function Sol

end module ode_mod

program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use farkode_mod           ! Fortran interface to the ARKode module
  use farkode_arkstep_mod   ! Fortran interface to the ARKStep module
  use fnvector_complex_mod  ! Custom complex N_Vector
  use ode_mod               ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  integer(c_int) :: ierr, iout
  real(c_double) :: tcur(1), tout, dTout, yerr, yerrI, yerr2

  type(N_Vector), pointer :: sunvec_y    ! sundials vector
  type(c_ptr)             :: arkode_mem  ! ARKODE memory

  ! solution vector
  type(FVec), pointer :: y

  !======= Internals ============

  ! initial problem output
  print *, "  "
  print *, "Analytical ODE test problem:"
  print '(2(a,f5.2),a)', "    lambda = (", real(lambda), " , ", imag(lambda), " ) "
  print '(2(a,es8.1))', "    reltol = ",reltol,",  abstol = ",abstol

  ! initialize SUNDIALS solution vector
  sunvec_y => FN_VNew_Complex(neq)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  y => FN_VGetFVec(sunvec_y)

  ! Set initial conditions into y
  y%data(1) = Sol(T0)

  ! create ARKStep memory
  arkode_mem = FARKStepCreate(c_funloc(Rhs), c_null_funptr, T0, sunvec_y)
  if (.not. c_associated(arkode_mem)) then
     print *,'ERROR: arkode_mem = NULL'
     stop 1
  end if

  ! main time-stepping loop: calls FARKStepEvolve to perform the integration, then
  ! prints results.  Stops when the final time has been reached
  tcur(1) = T0
  dTout = (Tf-T0)/Nt
  tout = T0+dTout
  yerrI = 0.d0
  yerr2 = 0.d0
  print *, " "
  print *, "      t     real(u)    imag(u)    error"
  print *, "   -------------------------------------------"
  print '(5x,f4.1,2(2x,es9.2),2x,es8.1)', tcur(1), real(y%data(1)), imag(y%data(1)), 0.d0
  do iout = 1,Nt

     ! call integrator
     ierr = FARKStepEvolve(arkode_mem, tout, sunvec_y, tcur, ARK_NORMAL)
     if (ierr /= 0) then
        write(*,*) 'Error in FARKStepEvolve, ierr = ', ierr, '; halting'
        stop 1
     endif

     ! compute/accumulate solution error
     yerr  = abs( y%data(1) - Sol(tcur(1)) )
     yerrI = max(yerrI, yerr)
     yerr2 = yerr2 + yerr**2

     ! print solution statistics
     print '(5x,f4.1,2(2x,es9.2),2x,es8.1)', tcur(1), real(y%data(1)), imag(y%data(1)), yerr

     ! update output time
     tout = min(tout + dTout, Tf)

  end do
  yerr2 = dsqrt( yerr2 / Nt )
  print *, "   -------------------------------------------"

  ! diagnostics output
  call ARKStepStats(arkode_mem)
  print '(2(a,es9.2))', "    Error: max = ", yerrI, ", rms = ", yerr2
  print *, ' '

  ! clean up
  call FARKStepFree(arkode_mem)
  call FN_VDestroy(sunvec_y)

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

  integer(c_int)  :: ierr         ! error flag
  integer(c_long) :: nsteps(1)    ! num steps
  integer(c_long) :: nst_a(1)     ! num steps attempted
  integer(c_long) :: nfe(1)       ! num explicit function evals
  integer(c_long) :: nfi(1)       ! num implicit function evals
  integer(c_long) :: netfails(1)  ! num error test fails

  !======= Internals ============

  ierr = FARKStepGetNumSteps(arkode_mem, nsteps)
  ierr = FARKStepGetNumStepAttempts(arkode_mem, nst_a)
  ierr = FARKStepGetNumRhsEvals(arkode_mem, nfe, nfi)
  ierr = FARKStepGetNumErrTestFails(arkode_mem, netfails)

  print *, ' '
  print *, 'Final Solver Statistics:'
  print '(4x,2(A,i4),A)' ,'Internal solver steps = ',nsteps(1),', (attempted = ',nst_a(1),')'
  print '(4x,A,i5)'   ,'Total RHS evals = ',nfe(1)
  print '(4x,A,i5)'      ,'Total number of error test failures =',netfails(1)

  return

end subroutine ARKStepStats
