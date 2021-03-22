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
! The following is a simple example problem for IDA, due to Robertson,
! is from chemical kinetics, and consists of the following three
! equations:
!
!      dy1/dt = -.04*y1 + 1.e4*y2*y3
!      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!         0   = y1 + y2 + y3 - 1
!
! on the interval from t = 0.0 to t = 4.e10, with initial
! conditions: y1 = 1, y2 = y3 = 0.
!
! While integrating the system, we also use the rootfinding
! feature to find the points at which y1 = 1e-4 or at which
! y3 = 0.01.
!
! The problem is solved with IDA using the DENSE linear
! solver, with a user-supplied Jacobian. Output is printed at
! t = .4, 4, 40, ..., 4e10.
! ------------------------------------------------------------------

module dae_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none

  integer(c_long), parameter :: neq = 3
  integer(c_long), parameter :: nout = 12

contains

  ! ----------------------------------------------------------------
  ! resrob: The DAE residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function resrob(tres, sunvec_y, sunvec_yp, sunvec_r, user_data) &
       result(ierr) bind(C,name='resrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tres      ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_yp ! derivative N_Vector
    type(N_Vector)        :: sunvec_r  ! residual N_Vector
    type(c_ptr),    value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yval(:)
    real(c_double), pointer :: ypval(:)
    real(c_double), pointer :: rval(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yval  => FN_VGetArrayPointer(sunvec_y)
    ypval => FN_VGetArrayPointer(sunvec_yp)
    rval  => FN_VGetArrayPointer(sunvec_r)

    ! fill residual vector
    rval(1)  = -0.04d0*yval(1) + 1.0d4*yval(2)*yval(3)
    rval(2)  = -rval(1) - 3.0d7*yval(2)**2 - ypval(2)
    rval(1)  = rval(1) - ypval(1)
    rval(3)  = yval(1) + yval(2) + yval(3) - 1.d0

    ! return success
    ierr = 0
    return

  end function resrob

  ! ----------------------------------------------------------------
  ! grob: The root function routine
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function grob(t, sunvec_y, sunvec_yp, gout, user_data) &
       result(ierr) bind(C,name='grob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_yp ! derivative N_Vector
    real(c_double)        :: gout(2)   ! root function values
    type(c_ptr),    value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yval(:)

    !======= Internals ============

    ! get data array from SUNDIALS vector
    yval => FN_VGetArrayPointer(sunvec_y)

    ! fill root vector
    gout(1) = yval(1) - 0.0001d0
    gout(2) = yval(3) - 0.01d0

    ! return success
    ierr = 0
    return

  end function grob

  ! ----------------------------------------------------------------
  ! jacrob: The DAE Jacobian function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function jacrob(t, cj, sunvec_y, sunvec_yp, sunvec_r, &
       sunmat_J, user_data, sunvec_t1, sunvec_t2, sunvec_t3) &
       result(ierr) bind(C,name='jacrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    use fsundials_matrix_mod
    use fnvector_serial_mod
    use fsunmatrix_dense_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: t         ! current time
    real(c_double), value :: cj        ! step size scaling factor
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_yp ! derivative N_Vector
    type(N_Vector)        :: sunvec_r  ! residual N_Vector
    type(SUNMatrix)       :: sunmat_J  ! Jacobian SUNMatrix
    type(c_ptr),    value :: user_data ! user-defined data
    type(N_Vector)        :: sunvec_t1 ! temporary N_Vectors
    type(N_Vector)        :: sunvec_t2
    type(N_Vector)        :: sunvec_t3

    ! pointers to data in SUNDIALS vector and matrix
    real(c_double), pointer :: yval(:)
    real(c_double), pointer :: J(:,:)


    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    j(1:3, 1:3) => FSUNDenseMatrix_Data(sunmat_J)

    ! fill Jacobian entries
    J(1,1) = -0.04d0 - cj
    J(2,1) = 0.04d0
    J(3,1) = 1.d0
    J(1,2) = 1.d4*yval(3)
    J(2,2) = -1.d4*yval(3) - 6.0d7*yval(2) - cj
    J(3,2) = 1.d0
    J(1,3) = 1.d4*yval(2)
    J(2,3) = -1.d4*yval(2)
    J(3,3) = 1.d0

    ! return success
    ierr = 0
    return

  end function jacrob

  ! ----------------------------------------------------------------
  ! check_ans: checks the solution error
  ! ----------------------------------------------------------------
  integer(c_int) function check_ans(y, t, rtol, atol) result(passfail)

    !======= Inclusions ===========
    use iso_c_binding

    !======= Declarations =========
    implicit none
    real(c_double) :: y(neq), atol(neq), t, rtol
    real(c_double) :: ref(neq), ewt(neq), err

    !======= Internals ============

    ! set the reference solution data
    ref(1) = 5.2083474251394888d-8
    ref(2) = 2.0833390772616859d-13
    ref(3) = 9.9999994791631752d-1

    ! compute the error weight vector, loosen atol by 10x
    ewt = 1.d0/(rtol*dabs(ref) + 10.d0*atol)

    ! compute the solution error
    ref = y-ref
    err = dsqrt( dot_product(ewt*ref,ewt*ref)/3 )

    ! is the solution within the tolerances (pass=0 or fail=1)?
    passfail = 0
    if (err .ge. 1.d0) then
       passfail = 1
       print *, " "
       print *, "SUNDIALS_WARNING: check_ans error=", err
       print *, " "
    end if

    return

  end function check_ans

end module dae_mod
! ------------------------------------------------------------------


program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  use fida_mod                      ! Fortran interface to IDA
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver
  use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver
  use fsundials_matrix_mod          ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod         ! Fortran interface to generic N_Vector
  use fsundials_linearsolver_mod    ! Fortran interface to generic SUNLinearSolver
  use fsundials_nonlinearsolver_mod ! Fortran interface to generic SUNNonlinearSolver
  use dae_mod                       ! ODE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: rtol, t0, tout1, tout, tret(1)
  integer(c_int) :: iout, retval, retvalr, nrtfn, rootsfound(2)

  type(N_Vector),           pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector),           pointer :: sunvec_yp     ! sundials derivative vector
  type(N_Vector),           pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix),          pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem       ! IDA memory

  ! solution and tolerance vectors, neq is set in the dae_mod module
  real(c_double) :: yval(neq), ypval(neq), avtol(neq)

  !======= Internals ============

  ! initialize solution vectors and tolerances
  yval(1) = 1.d0
  yval(2) = 0.d0
  yval(3) = 0.d0

  ypval(1) = -0.040
  ypval(2) = 0.04d0
  ypval(3) = 0.d0

  rtol = 1.d-4

  avtol(1) = 1.d-8
  avtol(2) = 1.d-6
  avtol(3) = 1.d-6

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(neq, yval)
  if (.not. associated(sunvec_y)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_yp => FN_VMake_Serial(neq, ypval)
  if (.not. associated(sunvec_yp)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sunvec_av => FN_VMake_Serial(neq, avtol)
  if (.not. associated(sunvec_av)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  ! set integration limits
  t0 = 0.d0
  tout1 = 0.4d0

  call PrintHeader(rtol, avtol, yval)

  ! Call FIDACreate and FIDAInit to initialize IDA memory
  ida_mem = FIDACreate()
  if (.not. c_associated(ida_mem)) then
     print *, 'ERROR: ida_mem = NULL'
     stop 1
  end if

  retval = FIDAInit(ida_mem, c_funloc(resrob), t0, sunvec_y, sunvec_yp)
  if (retval /= 0) then
     print *, 'Error in FIDAInit, retval = ', retval, '; halting'
     stop 1
  end if

  ! Call FIDASVtolerances to set tolerances
  retval = FIDASVtolerances(ida_mem, rtol, sunvec_av)
  if (retval /= 0) then
     print *, 'Error in FIDASVtolerances, retval = ', retval, '; halting'
     stop 1
  end if

  ! Call FIDARootInit to specify the root function grob with 2 components
  nrtfn = 2
  retval = FIDARootInit(ida_mem, nrtfn, c_funloc(grob))
  if (retval /= 0) then
     print *, 'Error in FIDARootInit, retval = ', retval, '; halting'
     stop 1
  end if

  ! Create dense SUNMatrix for use in linear solves
  sunmat_A => FSUNDenseMatrix(neq, neq)
  if (.not. associated(sunmat_A)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if

  ! Create dense SUNLinearSolver object
  sunlinsol_LS => FSUNDenseLinearSolver(sunvec_y, sunmat_A)
  if (.not. associated(sunlinsol_LS)) then
     print *, 'ERROR: sunlinsol = NULL'
     stop 1
  end if

  ! Attach the matrix and linear solver
  retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A);
  if (retval /= 0) then
     print *, 'Error in FIDASetLinearSolver, retval = ', retval, '; halting'
     stop 1
  end if

  ! Set the user-supplied Jacobian routine
  retval = FIDASetJacFn(ida_mem, c_funloc(jacrob))
  if (retval /= 0) then
     print *, 'Error in FIDASetJacFn, retval = ', retval, '; halting'
     stop 1
  end if

  ! Create Newton SUNNonlinearSolver object. IDA uses a
  ! Newton SUNNonlinearSolver by default, so it is not necessary
  ! to create it and attach it. It is done in this example code
  ! solely for demonstration purposes.
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y)
  if (.not. associated(sunnonlin_NLS)) then
     print *, 'ERROR: sunnonlinsol = NULL'
     stop 1
  end if

  ! Attach the nonlinear solver
  retval = FIDASetNonlinearSolver(ida_mem, sunnonlin_NLS)
  if (retval /= 0) then
     print *, 'Error in FIDASetNonlinearSolver, retval = ', retval, '; halting'
     stop 1
  end if

  ! In loop, call IDASolve, print results, and test for error.
  ! Break out of loop when NOUT preset output times have been reached.

  iout = 0
  tout = tout1
  do

     retval = FIDASolve(ida_mem, tout, tret, sunvec_y, sunvec_yp, IDA_NORMAL)
     if (retval < 0) then
        print *, 'Error in FIDASolve, retval = ', retval, '; halting'
        stop 1
     endif

     call PrintOutput(ida_mem, tret(1), yval)

     if (retval .eq. IDA_ROOT_RETURN) then
        retvalr = FIDAGetRootInfo(ida_mem, rootsfound)
        if (retvalr < 0) then
           print *, 'Error in FIDAGetRootInfo, retval = ', retval, '; halting'
           stop 1
        endif
        print '(a,2(i2,2x))', "    rootsfound[] = ", rootsfound(1), rootsfound(2)
     end if

     if (retval .eq. IDA_SUCCESS) then
        iout = iout+1
        tout = tout*10.d0
     end if

     if (iout .eq. NOUT) exit

  end do

  call PrintFinalStats(ida_mem)

  retval = check_ans(yval, tret(1), rtol, avtol)

  ! free memory
  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  call FN_VDestroy(sunvec_yp)

end program main


! ----------------------------------------------------------------
! PrintHeader: prints first lines of output (problem description)
! ----------------------------------------------------------------
subroutine PrintHeader(rtol, avtol, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use dae_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double) :: rtol
  real(c_double) :: avtol(neq)
  real(c_double) :: y(neq)

  !======= Internals ============

  print *, " "
  print *, "idaRoberts_dns_f2003: Robertson kinetics DAE serial example problem for IDA"
  print *, "         Three equation chemical kinetics problem."
  print *, " "
  print *, "Linear solver: DENSE, with user-supplied Jacobian."
  print '(a,f6.4,a,3(es6.0,1x))', "Tolerance parameters:  rtol = ",rtol,"   atol = ", avtol
  print '(a,3(f5.2,1x),a)', "Initial conditions y0 = (",y,")"
  print *, "Constraints and id not used."
  print *, " "
  print *, "-----------------------------------------------------------------------"
  print *, "  t             y1           y2           y3     | nst  k      h"
  print *, "-----------------------------------------------------------------------"

  return
end subroutine PrintHeader


! ----------------------------------------------------------------
! PrintOutput
! ----------------------------------------------------------------
subroutine PrintOutput(ida_mem, t, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fida_mod
  use dae_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  type(c_ptr)    :: ida_mem
  real(c_double) :: t, y(neq)

  ! internal variables
  integer(c_int)  :: retval, kused(1)
  integer(c_long) :: nst(1)
  real(c_double)  :: hused(1)

  !======= Internals ============

  retval = FIDAGetLastOrder(ida_mem, kused)
  if (retval /= 0) then
     print *, 'Error in FIDAGetLastOrder, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumSteps(ida_mem, nst)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetLastStep(ida_mem, hused)
  if (retval /= 0) then
     print *, 'Error in FIDAGetLastStep, retval = ', retval, '; halting'
     stop 1
  end if

  print '(es10.4,1x,3(es12.4,1x),a,i3,2x,i1,1x,es12.4)', &
       t, y(1), y(2), y(3), "| ", nst, kused(1), hused(1)

end subroutine PrintOutput


! ----------------------------------------------------------------
! PrintFinalStats
!
! Print KINSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(ida_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use fida_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: ida_mem
  integer(c_int)  :: retval
  integer(c_long) :: nst(1), nni(1), nje(1), nre(1), nreLS(1), nge(1), ncfn(1), netf(1)

  !======= Internals ============

  retval = FIDAGetNumSteps(ida_mem, nst)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumResEvals(ida_mem, nre)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumResEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumJacEvals(ida_mem, nje)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumJacEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumNonlinSolvIters(ida_mem, nni)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumNonlinSolvIters, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumErrTestFails(ida_mem, netf)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumErrTestFails, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumNonlinSolvConvFails(ida_mem, ncfn)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumLinResEvals(ida_mem, nreLS)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumLinResEvals, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FIDAGetNumGEvals(ida_mem, nge)
  if (retval /= 0) then
     print *, 'Error in FIDAGetNumGEvals, retval = ', retval, '; halting'
     stop 1
  end if

  print *, " "
  print *, "Final Run Statistics: "
  print *, "Number of steps                    = ", nst
  print *, "Number of residual evaluations     = ", nre+nreLS
  print *, "Number of Jacobian evaluations     = ", nje
  print *, "Number of nonlinear iterations     = ", nni
  print *, "Number of error test failures      = ", netf
  print *, "Number of nonlinear conv. failures = ", ncfn
  print *, "Number of root fn. evaluations     = ", nge

  return

end subroutine PrintFinalStats
