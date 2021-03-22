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
! Program to test custom fnvector_fortran_mod implementation
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Utility module for error-checking
! ------------------------------------------------------------------
module fnvector_test_mod
  use, intrinsic :: iso_c_binding
  use fnvector_fortran_mod
  implicit none

contains
  integer(c_int) function check_ans(val, tol, Nvar, N, sunvec_x) result(failure)

    implicit none
    real(c_double), value :: val, tol
    integer(c_long), value :: Nvar, N
    Type(N_Vector) :: sunvec_x
    Type(FVec), pointer :: x
    integer(c_long) :: i, j

    x => FN_VGetFVec(sunvec_x)
    failure = 0
    do j = 1,N
       do i = 1,Nvar
          if (dabs(x%data(i,j) - val) > tol)  failure = 1
       end do
    end do

  end function check_ans
end module fnvector_test_mod

! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fnvector_fortran_mod
  use fnvector_test_mod

  !======= Declarations =========
  implicit none

  ! local variables
  integer(c_int)             :: fails
  integer(c_long)            :: i, j, loc
  integer(c_long), parameter :: N = 1000
  integer(c_long), parameter :: Nvar = 10
  type(N_Vector), pointer    :: sU, sV, sW, sX, sY, sZ
  type(FVec), pointer        :: U, V, W, X, Y, Z
  real(c_double)             :: Udata(Nvar,N)
  real(c_double)             :: fac
  logical                    :: failure


  !======= Internals ============

  ! initialize failure total
  fails = 0

  ! create new vectors, using New, Make and Clone routines
  sU => FN_VMake_Fortran(Nvar, N, Udata)
  if (.not. associated(sU)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  U => FN_VGetFVec(sU)

  sV => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sV)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  V => FN_VGetFVec(sV)

  sW => FN_VNew_Fortran(Nvar, N)
  if (.not. associated(sW)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  W => FN_VGetFVec(sW)

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

  call c_f_pointer(FN_VClone_Fortran(sU), sZ)
  if (.not. associated(sZ)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  Z => FN_VGetFVec(sZ)


  ! check vector ID
  if (FN_VGetVectorID(sU) /= SUNDIALS_NVEC_CUSTOM) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VGetVectorID'
     print *, '    Unrecognized vector type', FN_VGetVectorID(sU)
  else
     print *, 'PASSED test -- FN_VGetVectorID'
  end if


  ! check vector length
  if (FN_VGetLength(sV) /= (N*Nvar)) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VGetLength'
     print *, '    ', FN_VGetLength(sV), ' /= ', N*Nvar
  else
     print *, 'PASSED test -- FN_VGetLength'
  end if

  ! test FN_VConst
  Udata = 0.d0
  call FN_VConst(1.d0, sU)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sU) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VConst'
  else
     print *, 'PASSED test -- FN_VConst'
  end if

  ! test FN_VLinearSum
  call FN_VConst(1.d0, sX)
  call FN_VConst(-2.d0, sY)
  call FN_VLinearSum(1.d0, sX, 1.d0, sY, sY)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sY) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 1a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 1a'
  end if

  call FN_VConst(1.d0, sX)
  call FN_VConst(2.d0, sY)
  call FN_VLinearSum(-1.d0, sX, 1.d0, sY, sY)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sY) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 1b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 1b'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(-2.d0, sY)
  call FN_VLinearSum(0.5d0, sX, 1.d0, sY, sY)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sY) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 1c'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 1c'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(-1.d0, sY)
  call FN_VLinearSum(1.d0, sX, 1.d0, sY, sX)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 2a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 2a'
  end if

  call FN_VConst(1.d0, sX)
  call FN_VConst(2.d0, sY)
  call FN_VLinearSum(1.d0, sX, -1.d0, sY, sX)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 2b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 2b'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(-0.5d0, sY)
  call FN_VLinearSum(1.d0, sX, 2.d0, sY, sX)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 2c'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 2c'
  end if

  call FN_VConst(-2.d0, sX)
  call FN_VConst(1.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(1.d0, sX, 1.d0, sY, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 3'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 3'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(1.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(1.d0, sX, -1.d0, sY, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 4a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 4a'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(1.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(-1.d0, sX, 1.d0, sY, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 4b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 4b'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(-0.5d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(1.d0, sX, 2.d0, sY, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 5a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 5a'
  end if

  call FN_VConst(0.5d0, sX)
  call FN_VConst(-2.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(2.d0, sX, 1.d0, sY, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 5b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 5b'
  end if

  call FN_VConst(-2.d0, sX)
  call FN_VConst(-0.5d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(-1.d0, sX, 2.d0, sY, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 6a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 6a'
  end if

  call FN_VConst(0.5d0, sX)
  call FN_VConst(2.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(2.d0, sX, -1.d0, sY, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 6b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 6b'
  end if

  call FN_VConst(1.d0, sX)
  call FN_VConst(-0.5d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(2.d0, sX, 2.d0, sY, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 7'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 7'
  end if

  call FN_VConst(0.5d0, sX)
  call FN_VConst(1.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(2.d0, sX, -2.d0, sY, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 8'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 8'
  end if

  call FN_VConst(1.d0, sX)
  call FN_VConst(-2.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VLinearSum(2.d0, sX, 0.5d0, sY, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 9'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 9'
  end if

  ! test FN_VProd
  call FN_VConst(2.d0, sX)
  call FN_VConst(-0.5d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VProd(sX, sY, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VProd'
  else
     print *, 'PASSED test -- FN_VProd'
  end if

  ! test FN_VDiv
  call FN_VConst(1.d0, sX)
  call FN_VConst(2.d0, sY)
  call FN_VConst(0.d0, sZ)
  call FN_VDiv(sX, sY, sZ)
  if (check_ans(0.5d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VDiv'
  else
     print *, 'PASSED test -- FN_VDiv'
  end if

  ! test FN_VScale
  call FN_VConst(0.5d0, sX)
  call FN_VScale(2.d0, sX, sX)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 1'
  else
     print *, 'PASSED test -- FN_VScale Case 1'
  end if

  call FN_VConst(-1.d0, sX)
  call FN_VConst(0.d0, sZ)
  call FN_VScale(1.d0, sX, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 2'
  else
     print *, 'PASSED test -- FN_VScale Case 2'
  end if

  call FN_VConst(-1.d0, sX)
  call FN_VConst(0.d0, sZ)
  call FN_VScale(-1.d0, sX, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 3'
  else
     print *, 'PASSED test -- FN_VScale Case 3'
  end if

  call FN_VConst(-0.5d0, sX)
  call FN_VConst(0.d0, sZ)
  call FN_VScale(2.d0, sX, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 4'
  else
     print *, 'PASSED test -- FN_VScale Case 4'
  end if

  ! test FN_VAbs
  call FN_VConst(-1.d0, sX)
  call FN_VConst(0.d0, sZ)
  call FN_VAbs(sX, sZ)
  if (check_ans(1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VAbs'
  else
     print *, 'PASSED test -- FN_VAbs'
  end if

  ! test FN_VInv
  call FN_VConst(2.d0, sX)
  call FN_VConst(0.d0, sZ)
  call FN_VInv(sX, sZ)
  if (check_ans(0.5d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInv'
  else
     print *, 'PASSED test -- FN_VInv'
  end if

  ! test FN_VAddConst
  call FN_VConst(1.d0, sX)
  call FN_VConst(0.d0, sZ)
  call FN_VAddConst(sX, -2.d0, sZ)
  if (check_ans(-1.d0, 1.d-14, Nvar, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VAddConst'
  else
     print *, 'PASSED test -- FN_VAddConst'
  end if

  ! test FN_VDotProd
  call FN_VConst(2.d0, sX)
  call FN_VConst(0.5d0, sY)
  if (dabs(FN_VDotProd(sX,sY) - (N*Nvar)) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VDotProd (',FN_VDotProd(sX,sY),' /= ',N*Nvar,')'
  else
     print *, 'PASSED test -- FN_VDotProd'
  end if

  ! test FN_VMaxNorm
  call FN_VConst(-0.5d0, sX)
  X%data(Nvar,N) = -2.d0
  if (dabs(FN_VMaxNorm(sX) - 2.d0) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VMaxNorm (',FN_VMaxNorm(sX),' /= 2.d0)'
  else
     print *, 'PASSED test -- FN_VMaxNorm'
  end if

  ! test FN_VWrmsNorm
  call FN_VConst(-0.5d0, sX)
  call FN_VConst(0.5d0, sY)
  if (dabs(FN_VWrmsNorm(sX,sY) - 0.25d0) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWrmsNorm (',FN_VWrmsNorm(sX,sY),' /= 0.25d0)'
  else
     print *, 'PASSED test -- FN_VWrmsNorm'
  end if

  ! test FN_VWrmsNormMask
  call FN_VConst(-0.5d0, sX)
  call FN_VConst(0.5d0, sY)
  call FN_VConst(1.d0, sZ)
  Z%data(Nvar,N) = 0.d0
  fac = dsqrt(1.d0*(N*Nvar - 1)/(N*Nvar))*0.25d0
  if (dabs(FN_VWrmsNormMask(sX,sY,sZ) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWrmsNormMask (',FN_VWrmsNormMask(sX,sY,sZ),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWrmsNormMask'
  end if

  ! test FN_VMin
  call FN_VConst(2.d0, sX)
  X%data(Nvar,N) = -2.d0
  if (dabs(FN_VMin(sX) + 2.d0) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VMin (',FN_VMin(sX),' /= -2.d0)'
  else
     print *, 'PASSED test -- FN_VMin'
  end if

  ! test FN_VWL2Norm
  call FN_VConst(-0.5d0, sX)
  call FN_VConst(0.5d0, sY)
  fac = dsqrt(1.d0*N*Nvar)*0.25d0
  if (dabs(FN_VWL2Norm(sX,sY) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWL2Norm (',FN_VWL2Norm(sX,sY),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWL2Norm'
  end if

  ! test FN_VL1Norm
  call FN_VConst(-1.d0, sX)
  fac = 1.d0*N*Nvar
  if (dabs(FN_VL1Norm(sX) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VL1Norm (',FN_VL1Norm(sX),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VL1Norm'
  end if

  ! test FN_VCompare
  call FN_VConst(-1.d0, sZ)
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 3)
        if (loc == 0)  X%data(i,j) = 0.d0
        if (loc == 1)  X%data(i,j) = -1.d0
        if (loc == 2)  X%data(i,j) = -2.d0
     end do
  end do
  call FN_VCompare(1.d0, sX, sZ)
  failure = .false.
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 3)
        if ((loc == 0) .and. (Z%data(i,j) /= 0.d0))  failure = .true.
        if ((loc == 1) .and. (Z%data(i,j) /= 1.d0))  failure = .true.
        if ((loc == 2) .and. (Z%data(i,j) /= 1.d0))  failure = .true.
     end do
  end do
  if (failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VCompare'
  else
     print *, 'PASSED test -- FN_VCompare'
  end if

  ! test FN_VInvTest
  call FN_VConst(0.5d0, sX)
  call FN_VConst(0.d0, sZ)
  failure = (FN_VInvTest(sX, sZ) == 0)
  if ((check_ans(2.d0, 1.d-14, Nvar, N, sZ) /= 0) .or. failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInvTest Case 1'
  else
     print *, 'PASSED test -- FN_VInvTest Case 1'
  end if

  failure = .false.
  call FN_VConst(0.d0, sZ)
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 2)
        if (loc == 0)  X%data(i,j) = 0.d0
        if (loc == 1)  X%data(i,j) = 0.5d0
     end do
  end do
  if (FN_VInvTest(sX, sZ) == 1)  failure = .true.
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 2)
        if ((loc == 0) .and. (Z%data(i,j) /= 0.d0))  failure = .true.
        if ((loc == 1) .and. (Z%data(i,j) /= 2.d0))  failure = .true.
     end do
  end do
  if (failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInvTest Case 2'
  else
     print *, 'PASSED test -- FN_VInvTest Case 2'
  end if

  ! test FN_VConstrMask
  call FN_VConst(-1.d0, sZ)
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 7)
        if (loc == 0) then  ! y = -2, test for < 0
           Y%data(i,j) = -2.d0
           X%data(i,j) = -2.d0
        end if
        if (loc == 1) then ! y = -1, test for <= 0
           Y%data(i,j) = -1.d0
           X%data(i,j) = -1.d0
        end if
        if (loc == 2) then ! y = -1, test for == 0
           Y%data(i,j) = -1.d0
           X%data(i,j) = 0.d0
        end if
        if (loc == 3) then ! y = 0, no test
           Y%data(i,j) = 0.d0
           X%data(i,j) = 0.5d0
        end if
        if (loc == 4) then ! y = 1, test for == 0
           Y%data(i,j) = 1.d0
           X%data(i,j) = 0.d0
        end if
        if (loc == 5) then ! y = 1, test for >= 0
           Y%data(i,j) = 1.d0
           X%data(i,j) = 1.d0
        end if
        if (loc == 6) then ! y = 2, test for > 0
           Y%data(i,j) = 2.d0
           X%data(i,j) = 2.d0
        end if
     end do
  end do
  failure = .false.
  if (FN_VConstrMask(sY, sX, sZ) /= 1) then
     failure = .true.
  end if
  if ((check_ans(0.d0, 1.d-14, Nvar, N, sZ) /= 0) .or. failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VConstrMask Case 1'
  else
     print *, 'PASSED test -- FN_VConstrMask Case 1'
  end if

  call FN_VConst(-1.d0, sZ)
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 5)
        if (loc == 0) then  ! y = -2, test for < 0
           Y%data(i,j) = -2.d0
           X%data(i,j) = 2.d0
        end if
        if (loc == 1) then ! y = -1, test for <= 0
           Y%data(i,j) = -1.d0
           X%data(i,j) = 1.d0
        end if
        if (loc == 2) then ! y = 0, no test
           Y%data(i,j) = 0.d0
           X%data(i,j) = 0.5d0
        end if
        if (loc == 3) then ! y = 1, test for >= 0
           Y%data(i,j) = 1.d0
           X%data(i,j) = -1.d0
        end if
        if (loc == 4) then ! y = 2, test for > 0
           Y%data(i,j) = 2.d0
           X%data(i,j) = -2.d0
        end if
     end do
  end do
  failure = .false.
  if (FN_VConstrMask(sY, sX, sZ) /= 0) then
     failure = .true.
  end if
  do j = 1,N
     do i = 1,Nvar
        loc = mod((j-1)*Nvar + i - 1, 5)
        if (loc == 2) then
           if (Z%data(i,j) /= 0.d0) failure = .true.
        else
           if (Z%data(i,j) /= 1.d0) failure = .true.
        end if
     end do
  end do
  if (failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VConstrMask Case 2'
  else
     print *, 'PASSED test -- FN_VConstrMask Case 2'
  end if

  ! test FN_VMinQuotient
  call FN_VConst(2.d0, sX)
  call FN_VConst(2.d0, sY)
  X%data(Nvar,N) = 0.5d0
  fac = 0.25d0
  if (dabs(FN_VMinQuotient(sX,sY) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VMinQuotient Case 1'
  else
     print *, 'PASSED test -- FN_VMinQuotient Case 1'
  end if

  call FN_VConst(2.d0, sX)
  call FN_VConst(0.d0, sY)
  fac = 1.d307
  if (dabs(FN_VMinQuotient(sX,sY) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VMinQuotient Case 2'
  else
     print *, 'PASSED test -- FN_VMinQuotient Case 2'
  end if

  ! test FN_VWSqrSumLocal
  call FN_VConst(-1.d0, sX)
  call FN_VConst(0.5d0, sY)
  fac = 0.25d0*N*Nvar
  if (dabs(FN_VWSqrSumLocal(sX,sY) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWSqrSumLocal (',FN_VWSqrSumLocal(sX,sY),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWSqrSumLocal'
  end if


  ! test FN_VWSqrSumMaskLocal
  call FN_VConst(-1.d0, sX)
  call FN_VConst(0.5d0, sY)
  call FN_VConst(1.d0, sZ)
  Z%data(Nvar,N) = 0.d0
  fac = 0.25d0*(N*Nvar-1)
  if (dabs(FN_VWSqrSumMaskLocal(sX,sY,sZ) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWSqrSumMaskLocal (',FN_VWSqrSumMaskLocal(sX,sY,sZ),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWSqrSumMaskLocal'
  end if

  ! free vectors
  call FN_VDestroy(sU)
  call FN_VDestroy(sV)
  call FN_VDestroy(sW)
  call FN_VDestroy(sX)
  call FN_VDestroy(sY)
  call FN_VDestroy(sZ)

  ! print results
  if (fails > 0) then
     print '(a,i3,a)', 'FAIL: FNVector module failed ',fails,' tests'
     stop 1
  else
     print *, 'SUCCESS: FNVector module passed all tests'
  end if
  print *, '  '

end program main
