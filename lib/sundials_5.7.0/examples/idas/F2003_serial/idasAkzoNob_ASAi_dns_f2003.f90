! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! Acknowledgements: Based on idaAkzoNob_ASAi_dns.c by Radu Serban
! and Cosmin Petra @ LLNL.
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
! Adjoint sensitivity example problem
!
! This IVP is a stiff system of 6 non-linear DAEs of index 1. The
! problem originates from Akzo Nobel Central research in Arnhern,
! The Netherlands, and describes a chemical process in which 2
! species are mixed, while carbon dioxide is continuously added.
! See http://pitagora.dm.uniba.it/~testset/report/chemakzo.pdf
!
! FIDAS also computes the sensitivities with respect to initial
! conditions of the following quantity:
!   G = int_t0^t1 y1 dt
! The sensitivity of G is the solution of the adjoint system at t0.
! -----------------------------------------------------------------

module dae_mod
  use, intrinsic :: iso_c_binding
  implicit none

  ! problem parameters
  integer(c_long), parameter :: NEQ   = 6
  integer(c_long), parameter :: STEPS = 150
  real(c_double),  parameter :: T0    = 0.0d0
  real(c_double),  parameter :: TF    = 180.d0
  real(c_double),  parameter :: RTOL  = 1e-08
  real(c_double),  parameter :: ATOL  = 1e-10
  real(c_double),  parameter :: RTOLB = 1e-06
  real(c_double),  parameter :: ATOLB = 1e-08
  real(c_double),  parameter :: RTOLQ = 1e-10
  real(c_double),  parameter :: ATOLQ = 1e-12

  ! problem constants
  real(c_double) :: ZERO    = 0.0d0
  real(c_double) :: QUARTER = 0.25d0
  real(c_double) :: HALF    = 0.5d0
  real(c_double) :: ONE     = 1.0d0
  real(c_double) :: TWO     = 2.0d0
  real(c_double) :: FOUR    = 4.0d0
  real(c_double) :: EIGHT   = 8.0d0

  ! problem data
  real(c_double) :: k1, k2, k3, k4
  real(c_double) :: K, klA, Ks, pCO2, H

contains

  integer(c_int) function res(t, nv_yy, nv_yd, nv_resval, userdata) &
    result(retval) bind(C,name='res')
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    ! function arguments
    real(c_double), value :: t
    type(N_Vector)        :: nv_yy, nv_yd, nv_resval
    type(c_ptr)           :: userdata

    ! local variables
    real(c_double)          :: y1, y2, y3, y4, y5, y6
    real(c_double)          :: yd1, yd2, yd3, yd4, yd5
    real(c_double)          :: r1, r2, r3, r4, r5, Fin
    real(c_double), pointer :: yy(:), yd(:), resval(:)

    yy     => FN_VGetArrayPointer(nv_yy)
    yd     => FN_VGetArrayPointer(nv_yd)
    resval => FN_VGetArrayPointer(nv_resval)

    y1 = yy(1)
    y2 = yy(2)
    y3 = yy(3)
    y4 = yy(4)
    y5 = yy(5)
    y6 = yy(6)

    yd1 = yd(1)
    yd2 = yd(2)
    yd3 = yd(3)
    yd4 = yd(4)
    yd5 = yd(5)

    r1  = k1 * (y1**4) * sqrt(y2)
    r2  = k2 * y3 * y4
    r3  = k2/K * y1 * y5
    r4  = k3 * y1 * y4 * y4
    r5  = k4 * y6 * y6 * sqrt(y2)
    Fin = klA * ( pCO2/H - y2 )

    resval(1) = yd1 + TWO*r1 - r2 + r3 + r4
    resval(2) = yd2 + HALF*r1 + r4 + HALF*r5 - Fin
    resval(3) = yd3 - r1 + r2 - r3
    resval(4) = yd4 + r2 - r3 + TWO*r4
    resval(5) = yd5 - r2 + r3 - r5
    resval(6) = Ks*y1*y4 - y6

    retval = 0
    return
  end function res

  integer(c_int) function rhsQ(t, nv_yy, nv_yp, nv_qdot, userdata) &
    result(retval) bind(C,name='rhsQ')
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    ! function arguments
    real(c_double), value :: t
    type(N_Vector)        :: nv_yy, nv_yp, nv_qdot
    type(c_ptr)           :: userdata

    ! local variables
    real(c_double), pointer :: qdot(:), yy(:)

    qdot => FN_VGetArrayPointer(nv_qdot)
    yy   => FN_VGetArrayPointer(nv_yy)
    qdot(1) = yy(1)

    retval = 0
    return
  end function rhsQ

  integer(c_int) function resB(tt, nv_yy, nv_yp, nv_yyB, nv_ypB, nv_rrB, userdata) &
    result(retval) bind(C,name='resB')
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    ! function arguments
    real(c_double), value :: tt
    type(N_Vector)        :: nv_yy, nv_yp, nv_yyB, nv_ypB, nv_rrB
    type(c_ptr)           :: userdata

    ! local variables
    real(c_double) :: y1, y2, y3, y4, y5, y6
    real(c_double) :: yB1, yB2, yB3, yB4, yB5, yB6
    real(c_double) :: ypB1, ypB2, ypB3, ypB4, ypB5
    real(c_double) :: y2tohalf, y1to3, k2overK, tmp1, tmp2
    real(c_double), pointer :: yy(:), yyB(:), ypB(:), rrb(:)

    yy  => FN_VGetArrayPointer(nv_yy)
    yyB => FN_VGetArrayPointer(nv_yyB)
    ypB => FN_VGetArrayPointer(nv_ypB)
    rrB => FN_VGetArrayPointer(nv_rrB)

    y1 = yy(1)
    y2 = yy(2)
    y3 = yy(3)
    y4 = yy(4)
    y5 = yy(5)
    y6 = yy(6)

    yB1 = yyB(1)
    yB2 = yyB(2)
    yB3 = yyB(3)
    yB4 = yyB(4)
    yB5 = yyB(5)
    yB6 = yyB(6)

    ypB1 = ypB(1)
    ypB2 = ypB(2)
    ypB3 = ypB(3)
    ypB4 = ypB(4)
    ypB5 = ypB(5)

    y2tohalf = sqrt(y2)
    y1to3 = y1*y1*y1
    k2overK = k2/K

    tmp1 = k1* y1to3 * y2tohalf
    tmp2 = k3*y4*y4
    rrB(1) = 1 +  ypB1 - (EIGHT*tmp1 + k2overK*y5 + tmp2)*yB1 &
           - (TWO*tmp1+tmp2)*yB2 + (FOUR*tmp1+k2overK*y5)*yB3 &
           + k2overK*y5*(yB4-yB5) - TWO*tmp2*yB4 + Ks*y4*yB6

    tmp1 = k1 * y1*y1to3 * (y2tohalf/y2)
    tmp2 = k4 * y6*y6 * (y2tohalf/y2)
    rrB(2) = ypB2 - tmp1*yB1 - (QUARTER*tmp1 + QUARTER*tmp2 + klA)*yB2 &
           + HALF*tmp1*yB3 + HALF*tmp2*yB5

    rrB(3) = ypB3 + k2*y4*(yB1-yB3-yB4+yB5)

    tmp1 = k3*y1*y4
    tmp2 = k2*y3
    rrB(4) = ypB4 + (tmp2-TWO*tmp1)*yB1 - TWO*tmp1*yB2 - tmp2*yB3 &
           - (tmp2+FOUR*tmp1)*yB4 + tmp2*yB5 + Ks*y1*yB6

    rrB(5) = ypB5 - k2overK*y1*(yB1-yB3-yB4+yB5)

    rrB(6) = k4*y6*y2tohalf*(2*yB5-yB2) - yB6

    retval = 0
    return
  end function resB

  subroutine PrintOutput(nv_yB, nv_ypB)
    use, intrinsic :: iso_c_binding
    use fsundials_nvector_mod
    implicit none

    ! function arguments
    type(N_Vector) :: nv_yB, nv_ypB
    real(c_double), pointer :: yB(:), ypB(:)

    yB => FN_VGetArrayPointer(nv_yB)

    write(*,'(1x,A,es12.4)') "dG/dy0:         ", yB(1)
    write(*,'(1x,A,es12.4)') "                ", yB(2)
    write(*,'(1x,A,es12.4)') "                ", yB(3)
    write(*,'(1x,A,es12.4)') "                ", yB(4)
    write(*,'(1x,A,es12.4)') "                ", yB(5)
    write(*,'(1x,A,es12.4)') "                ", yB(6)
    write(*,*) "--------------------------------------------------------"
    write(*,*) ""

  end subroutine

end module dae_mod


! Main program
program main
  use, intrinsic :: iso_c_binding
  use fidas_mod                  ! Fortran interface to IDA
  use fnvector_serial_mod        ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
  use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
  use fsundials_matrix_mod       ! Fortran interface to generic SUNMatrix
  use fsundials_nvector_mod      ! Fortran interface to generic N_Vector
  use fsundials_linearsolver_mod ! Fortran interface to generic SUNLinearSolver
  use dae_mod                    ! DAE problem module
  implicit none

  ! Local variables
  type(c_ptr)                    :: mem, memB
  integer(c_int)                 :: ncheck(1), retval
  real(c_double)                 :: time(1)
  integer(c_long)                :: nst(1), nstB(1)
  integer(c_int)                 :: indexB(1)
  real(c_double),        pointer :: yy(:), yp(:), rr(:), q(:)
  type(N_Vector),        pointer :: nv_yy, nv_yp, nv_rr, nv_q
  real(c_double),        pointer :: yB(:), ypB(:)
  type(N_Vector),        pointer :: nv_yB, nv_ypB
  type(SUNMatrix),       pointer :: A, AB
  type(SUNLinearSolver), pointer :: LS, LSB

  ! Consistent IC for  y, y'.
  real(c_double) :: y01 = 0.444d0
  real(c_double) :: y02 = 0.00123d0
  real(c_double) :: y03 = 0.0d0
  real(c_double) :: y04 = 0.007d0
  real(c_double) :: y05 = 0.0d0

  write(*,*) ""
  write(*,*) "Adjoint Sensitivity Example for Akzo-Nobel Chemical Kinetics"
  write(*,*) "-------------------------------------------------------------"
  write(*,*) "Sensitivity of G = int_t0^tf (y1) dt with respect to IC."
  write(*,*) "-------------------------------------------------------------"
  write(*,*) ""

  ! Fill problem data with the appropriate values for coefficients.
  k1   = 18.7d0
  k2   = 0.58d0
  k3   = 0.09d0
  k4   = 0.42d0
  K    = 34.4d0
  klA  = 3.3d0
  Ks   = 115.83d0
  pCO2 = 0.9d0
  H    = 737.0d0

  ! Allocate N-vectors.
  nv_yy => FN_VNew_Serial(NEQ)
  nv_yp => FN_VNew_Serial(NEQ)

  ! Set IC
  yy    => FN_VGetArrayPointer(nv_yy)
  yy(1) = y01
  yy(2) = y02
  yy(3) = y03
  yy(4) = y04
  yy(5) = y05
  yy(6) = Ks * y01 * y04

  ! Get y' = - res(t0, y, 0)
  call FN_VConst(ZERO, nv_yp)

  nv_rr  => FN_VNew_Serial(NEQ)
  retval = res(T0, nv_yy, nv_yp, nv_rr, c_null_ptr)
  call FN_VScale(-ONE, nv_rr, nv_yp)
  call FN_VDestroy(nv_rr)

  ! Create and initialize q0 for quadratures.
  nv_q => FN_VNew_Serial(1_8)
  if (.not. associated(nv_q)) then
    write(*,*) 'ERROR: FN_VNew_Serial returned NULL'
    stop 1
  end if

  q    => FN_VGetArrayPointer(nv_q)
  if (.not. associated(q)) then
    write(*,*) 'ERROR: FN_VGetArrayPointer returned NULL'
    stop 1
  end if
  q(1) = ZERO

  ! Call FIDACreate and FIDAInit to initialize FIDA memory
  mem    = FIDACreate()
  if (.not. c_associated(mem)) then
    write(*,*) 'ERROR: FIDACreate returned NULL'
    stop 1
  end if

  retval = FIDAInit(mem, c_funloc(res), T0, nv_yy, nv_yp)
  call check_retval(retval, "FIDAInit")

  ! Set tolerances.
  retval = FIDASStolerances(mem, RTOL, ATOL)
  call check_retval(retval, "FIDASStolerances")

  ! Create dense SUNMatrix for use in linear solves
  A => FSUNDenseMatrix(NEQ, NEQ)
  if (.not. associated(A)) then
    write(*,*) 'ERROR: FSUNDenseMatrix returned NULL'
    stop 1
  end if

  ! Create dense SUNLinearSolver object
  LS => FSUNLinSol_Dense(nv_yy, A)
  if (.not. associated(LS)) then
    write(*,*) 'ERROR: FSUNLinSol_Dense returned NULL'
    stop 1
  end if

  ! Attach the matrix and linear solver
  retval = FIDASetLinearSolver(mem, LS, A)
  call check_retval(retval, "FIDASetLinearSolver")

  ! Initialize QUADRATURE(S).
  retval = FIDAQuadInit(mem, c_funloc(rhsQ), nv_q)
  call check_retval(retval, "FIDAQuadInit")

  ! Set tolerances and error control for quadratures.
  retval = FIDAQuadSStolerances(mem, RTOLQ, ATOLQ)
  call check_retval(retval, "FIDAQuadSStolerances")

  retval = FIDASetQuadErrCon(mem, 1)
  call check_retval(retval, "FIDASetQuadErrCon")

  ! Prepare ADJOINT.
  retval = FIDAAdjInit(mem, STEPS, IDA_HERMITE)
  call check_retval(retval, "FIDAAdjInit")

  ! FORWARD run.
  write(*,'(1x,A)',advance='no') "Forward integration ... "
  retval = FIDASolveF(mem, TF, time, nv_yy, nv_yp, IDA_NORMAL, ncheck)
  call check_retval(retval, "FIDASolveF")
  retval = FIDAGetNumSteps(mem, nst)
  write(*,'(A,i6,A)') "done ( nst = ", nst, " )"
  retval = FIDAGetQuad(mem, time, nv_q)

  write(*,'(1x,A,F24.16)') "G:          ", q(1)
  write(*,*) "--------------------------------------------------------"
  write(*,*) ""

  ! BACKWARD run

  ! Initialize yB
  nv_yB => FN_VNew_Serial(NEQ)
  if (.not. associated(nv_yB)) then
    write(*,*) 'ERROR: FN_VNew_Serial returned NULL'
    stop 1
  end if
  call FN_VConst(ZERO, nv_yB)

  nv_ypB => FN_VNew_Serial(NEQ)
  if (.not. associated(nv_ypB)) then
    write(*,*) 'ERROR: FN_VNew_Serial returned NULL'
    stop 1
  end if
  call FN_VConst(ZERO, nv_ypB)
  ypB    => FN_VGetArrayPointer(nv_ypB)
  ypB(1) = -ONE

  retval = FIDACreateB(mem, indexB)
  call check_retval(retval, "FIDACreateB")

  retval = FIDAInitB(mem, indexB(1), c_funloc(resB), TF, nv_yB, nv_ypB)
  call check_retval(retval, "FIDAInitB")

  retval = FIDASStolerancesB(mem, indexB(1), RTOLB, ATOLB)
  call check_retval(retval, "FIDASStolerancesB")

  retval = FIDASetMaxNumStepsB(mem, indexB(1), 1000_8)
  call check_retval(retval, "FIDASetMaxNumStepsB")

  ! Create dense SUNMatrix for use in linear solves
  AB => FSUNDenseMatrix(NEQ, NEQ)
  if (.not. associated(AB)) then
    write(*,*) 'ERROR: FSUNDenseMatrix returned NULL'
    stop 1
  end if

  ! Create dense SUNLinearSolver object
  LSB => FSUNLinSol_Dense(nv_yB, AB)
  if (.not. associated(LSB)) then
    write(*,*) 'ERROR: FSUNLinSol_Dense returned NULL'
    stop 1
  end if

  ! Attach the matrix and linear solver
  retval = FIDASetLinearSolverB(mem, indexB(1), LSB, AB)
  call check_retval(retval, "FIDASetLinearSolverB")

  ! Do the backward integration
  write(*,'(1x,A)',advance='no') "Backward integration ... "
  retval = FIDASolveB(mem, T0, IDA_NORMAL)
  call check_retval(retval, "FIDASolveB")
  memB   = FIDAGetAdjIDABmem(mem, indexB(1))
  retval = FIDAGetNumSteps(memB, nstB)
  write(*,'(A,i6,A)') "done ( nst = ", nstB, " )"
  retval = FIDAGetB(mem, indexB(1), time, nv_yB, nv_ypB)

  ! Print the solution
  call PrintOutput(nv_yB, nv_ypB)

  ! Free memory
  call FIDAFree(mem)
  retval = FSUNLinSolFree(LS)
  retval = FSUNLinSolFree(LSB)
  call FSUNMatDestroy(A)
  call FSUNMatDestroy(AB)
  call FN_VDestroy(nv_yy)
  call FN_VDestroy(nv_yp)
  call FN_VDestroy(nv_yB)
  call FN_VDestroy(nv_ypB)
  call FN_VDestroy(nv_q)

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
