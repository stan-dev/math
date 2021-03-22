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
! Program to test custom fsunmatrix_fortran_mod implementation
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Utility module for error-checking
! ------------------------------------------------------------------
module fsunmatrix_test_mod
  use, intrinsic :: iso_c_binding
  use fsunmatrix_fortran_mod
  use fnvector_fortran_mod
  implicit none

contains
  ! ------------------------------------------------------------------
  integer(c_int) function check_matrix(sunmat_A, sunmat_B, tol, Nvar, N) result(failure)

    implicit none
    real(c_double), value :: tol
    integer(c_long), value :: Nvar, N
    Type(SUNMatrix) :: sunmat_A, sunmat_B
    Type(FMat), pointer :: A, B
    integer(c_long) :: i, j, k

    A => FSUNMatGetFMat(sunmat_A)
    B => FSUNMatGetFMat(sunmat_B)
    failure = 0
    do k = 1,N
       do j = 1,Nvar
          do i = 1,Nvar
             if (dabs(A%data(i,j,k) - B%data(i,j,k)) > tol)  failure = 1
          end do
       end do
    end do

  end function check_matrix

  ! ------------------------------------------------------------------
  integer(c_int) function check_matrix_entry(sunmat_A, val, tol, Nvar, N) result(failure)

    implicit none
    real(c_double), value :: tol, val
    integer(c_long), value :: Nvar, N
    Type(SUNMatrix) :: sunmat_A
    Type(FMat), pointer :: A
    integer(c_long) :: i, j, k

    A => FSUNMatGetFMat(sunmat_A)
    failure = 0
    do k = 1,N
       do j = 1,Nvar
          do i = 1,Nvar
             if (dabs(A%data(i,j,k) - val) > tol)  failure = 1
          end do
       end do
    end do

  end function check_matrix_entry

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
       print *, '    i      j        x        y       diff'
       print *, '  --------------------------------------------'
       do j = 1,N
          do i = 1,Nvar
             if (dabs(x%data(i,j) - y%data(i,j)) > tol) then
                print '(2x,2(i4,3x),3(es9.2,1x))', i, j, x%data(i,j), &
                     y%data(i,j), dabs(x%data(i,j) - y%data(i,j))
             end if
          end do
       end do
       print *, '  --------------------------------------------'
       print *, '  '
    end if

  end function check_vector

end module fsunmatrix_test_mod

! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fnvector_fortran_mod
  use fsunmatrix_test_mod

  !======= Declarations =========
  implicit none

  ! local variables
  integer(c_int)  :: fails, retval, i, j, k
  integer(c_long), parameter :: N = 1000
  integer(c_long), parameter :: Nvar = 50
  type(SUNMatrix), pointer :: sA, sB, sC, sD, sI
  type(FMat), pointer :: A, Eye
  type(N_Vector),  pointer :: sW, sX, sY, sZ
  type(FVec), pointer :: X, Y


  !======= Internals ============

  ! initialize failure total
  fails = 0

  ! create new matrices and vectors
  sW => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sW)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

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
  Y => FN_VGetFVec(sY)

  sZ => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sZ)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if

  sA => FSUNMatNew_Fortran(Nvar, N)
  if (.not. associated(sA)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if
  A => FSUNMatGetFMat(sA)

  sB => FSUNMatNew_Fortran(Nvar, N)
  if (.not. associated(sB)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if

  sC => FSUNMatNew_Fortran(Nvar, N)
  if (.not. associated(sC)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if

  sD => FSUNMatNew_Fortran(Nvar, N)
  if (.not. associated(sD)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if

  call c_f_pointer(FSUNMatClone_Fortran(sA), sI)
  if (.not. associated(sI)) then
     print *, 'ERROR: sunmat = NULL'
     stop 1
  end if
  Eye => FSUNMatGetFMat(sI)


  ! fill matrices and vectors
  X%data = 0.d0
  Y%data = 0.d0
  A%data = 0.d0
  Eye%data = 0.d0
  do k = 1,N
     do j = 1,Nvar
        do i = 1,Nvar
           A%data(i,j,k) = 1.d0*i*j/k
        end do
        Eye%data(j,j,k) = 1.d0
        x%data(j,k) = 1.d0*k/j
        y%data(j,k) = 1.d0*j*Nvar
     end do
  end do

  ! check matrix ID
  if (FSUNMatGetID(sA) /= SUNMATRIX_CUSTOM) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatGetID'
     print *, '    Unrecognized vector type', FSUNMatGetID(sA)
  else
     print *, 'PASSED test -- FSUNMatGetID'
  end if

  ! test SUNMatZero
  retval = FSUNMatZero(sB)
  if ( (check_matrix_entry(sB, 0.d0, 1.d-14, Nvar, N) /= 0) &
       .or. (retval /= SUNMAT_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatZero'
  else
     print *, 'PASSED test -- FSUNMatZero'
  end if

  ! test SUNMatCopy
  retval = FSUNMatCopy(sA, sB)
  if ( (check_matrix(sA, sB, 1.d-14, Nvar, N) /= 0) &
       .or. (retval /= SUNMAT_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatCopy'
  else
     print *, 'PASSED test -- FSUNMatCopy'
  end if

  ! test SUNMatScaleAdd
  retval = FSUNMatCopy(sA, sB)
  retval = FSUNMatScaleAdd(-1.d0, sB, sB)
  if ( (check_matrix_entry(sB, 0.d0, 1.d-14, Nvar, N) /= 0) &
       .or. (retval /= SUNMAT_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatScaleAdd case 1'
  else
     print *, 'PASSED test -- FSUNMatScaleAdd case 1'
  end if

  retval = FSUNMatCopy(sA, sD)
  retval = FSUNMatCopy(sI, sC)
  retval = FSUNMatScaleAdd(1.d0, sD, sI)
  if (retval == SUNMAT_SUCCESS)  retval = FSUNMatScaleAdd(1.d0, sC, sA)
  if ( (check_matrix(sD, sC, 1.d-14, Nvar, N) /= 0) &
       .or. (retval /= SUNMAT_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatScaleAdd case 2'
  else
     print *, 'PASSED test -- FSUNMatScaleAdd case 2'
  end if

  ! test SUNMatScaleAddI
  retval = FSUNMatCopy(sI, sB)
  retval = FSUNMatScaleAddI(-1.d0, sB)
  if ( (check_matrix_entry(sB, 0.d0, 1.d-14, Nvar, N) /= 0) &
       .or. (retval /= SUNMAT_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatScaleAddI'
  else
     print *, 'PASSED test -- FSUNMatScaleAddI'
  end if

  ! test SUNMatMatvec
  retval = FSUNMatCopy(sA, sB)
  retval = FSUNMatScaleAddI(3.d0, sB)
  retval = FSUNMatMatvec(sB, sX, sZ)
  call FN_VLinearSum(3.d0, sY, 1.d0, sX, sW)
  if ( (check_vector(sW, sZ, 1.d-15*Nvar*Nvar, Nvar, N) /= 0) &
       .or. (retval /= SUNMAT_SUCCESS) ) then
     fails = fails + 1
     print *, '>>> FAILED test -- FSUNMatMatvec'
  else
     print *, 'PASSED test -- FSUNMatMatvec'
  end if

  ! free matrices and vectors
  call FSUNMatDestroy(sA)
  call FSUNMatDestroy(sB)
  call FSUNMatDestroy(sC)
  call FSUNMatDestroy(sD)
  call FSUNMatDestroy(sI)
  call FN_VDestroy(sW)
  call FN_VDestroy(sX)
  call FN_VDestroy(sY)
  call FN_VDestroy(sZ)

  ! print results
  if (fails > 0) then
     print '(a,i3,a)', 'FAIL: FSUNMatrix module failed ',fails,' tests'
     stop 1
  else
     print *, 'SUCCESS: FSUNMatrix module passed all tests'
  end if
  print *, '  '

end program main
