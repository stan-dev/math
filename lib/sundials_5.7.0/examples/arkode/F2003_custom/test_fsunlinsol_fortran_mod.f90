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
! Program to test custom fsunlinsol_fortran_mod implementation
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Utility module for error-checking
! ------------------------------------------------------------------
module fsunlinsol_test_mod
  use, intrinsic :: iso_c_binding
  use fsunlinsol_fortran_mod
  use fsunmatrix_fortran_mod
  use fnvector_fortran_mod
  implicit none

contains
  ! ------------------------------------------------------------------
  integer(c_int) function check_vector(sunvec_x, sunvec_y, tol, Nvar, N) result(failure)

    implicit none
    real(c_double), value :: tol
    integer(c_long), value :: Nvar, N
    Type(N_Vector) :: sunvec_x, sunvec_y
    Type(FVec), pointer :: x, y
    integer(c_long) :: i, j

    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)
    failure = 0
    do j = 1,N
       do i = 1,Nvar
          if (dabs(x%data(i,j) - y%data(i,j)) > tol) then
             failure = 1
          end if
       end do
    end do

    if (failure == 1) then
       print *, '  '
       print *, 'check_vector failure, differences:'
       print *, '   blk    idx       x        y       diff'
       print *, '  --------------------------------------------'
       do j = 1,N
          do i = 1,Nvar
             if (dabs(x%data(i,j) - y%data(i,j)) > tol) then
                print '(2x,2(i4,3x),3(es9.2,1x))', j, i, x%data(i,j), &
                     y%data(i,j), dabs(x%data(i,j) - y%data(i,j))
             end if
          end do
       end do
       print *, '  --------------------------------------------'
       print *, '  '
    end if

  end function check_vector

end module fsunlinsol_test_mod

! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsunlinsol_test_mod

  !======= Declarations =========
  implicit none

  ! local variables
  integer(c_int)                   :: fails, retval, i, j, k
  integer(c_long),       parameter :: N = 1000
  integer(c_long),       parameter :: Nvar = 50
  type(SUNMatrix),       pointer   :: sA
  type(FMat),            pointer   :: A
  type(SUNLinearSolver), pointer   :: LS
  type(FLinSol),         pointer   :: S
  type(N_Vector),        pointer   :: sX, sY, sB
  type(FVec),            pointer   :: X, B


  !======= Internals ============

  ! initialize failure total
  fails = 0

  ! create new matrices and vectors
  sX => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sX)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  X => FN_VGetFVec(sX)

  sY => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sY)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sB => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sB)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  B => FN_VGetFVec(sB)

  sA => FSUNMatNew_Fortran(Nvar, N)
  if (.not. associated(sA)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if
  A => FSUNMatGetFMat(sA)


  ! fill A and X with uniformly-distributed random numbers in [0,1)
  call random_number(X%data)
  call random_number(A%data)

  ! update A to scale by 1/Nvar, and 1 to anti-diagonal of each diagonal block
  i = 0
  do k = 1,N
     A%data(:,:,k) = A%data(:,:,k)/Nvar
     do j = 1,Nvar
        A%data(Nvar-j+1,j,k) = A%data(Nvar-i+1,j,k) + 1.d0
     end do
  end do

  ! compute B = A*X
  retval = FSUNMatMatvec(sA, sX, sB)
  if (retval /= SUNMAT_SUCCESS) then
     print *, 'ERROR: FSUNMatMatvec fail'
     stop 1
  end if

  ! create custom linear solver
  LS => FSUNLinSolNew_Fortran(Nvar, N)
  if (.not. associated(LS)) then
     print *, 'ERROR: sunlinsol = NULL'
     stop 1
  end if
  S => FSUNLinSolGetFLinSol(LS)

  ! test SUNLinSolGetType
  if (FSUNLinSolGetType(LS) /= SUNLINEARSOLVER_DIRECT) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNLinSolGetType'
     print *, '    Unrecognized vector type', FSUNLinSolGetType(LS)
  else
     print *, 'PASSED test -- FSUNLinSolGetType'
  end if

  ! test SUNLinSolSetup
  retval = FSUNLinSolSetup(LS, sA)
  if (retval /= SUNLS_SUCCESS) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNLinSolSetup'
  else
     print *, 'PASSED test -- FSUNLinSolSetup'
  end if

  ! test SUNLinSolSolve
  call FN_VConst(0.d0, sY)
  retval = FSUNLinSolSolve(LS, sA, sY, sB, 1.d-9)
  if ( (check_vector(sX, sY, 1.d-15*Nvar*Nvar, Nvar, N) /= 0) &
       .or. (retval /= SUNLS_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNLinSolSolve'
  else
     print *, 'PASSED test -- FSUNLinSolSolve'
  end if

  ! free solver, matrix and vectors
  call FSUNMatDestroy(sA)
  call FN_VDestroy(sX)
  call FN_VDestroy(sY)
  call FN_VDestroy(sB)
  retval = FSUNLinSolFree(LS)
  if (retval /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNLinSolFree'
  else
     print *, 'PASSED test -- FSUNLinSolFree'
  end if

  ! print results
  if (fails > 0) then
     print '(a,i3,a)', 'FAIL: FSUNLinSol module failed ',fails,' tests'
     stop 1
  else
     print *, 'SUCCESS: FSUNLinSol module passed all tests'
  end if
  print *, '  '

end program main
