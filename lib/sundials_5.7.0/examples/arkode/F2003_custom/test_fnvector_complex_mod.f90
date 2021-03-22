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
! Program to test custom fnvector_complex_mod implementation
! ------------------------------------------------------------------

! ------------------------------------------------------------------
! Utility module for error-checking
! ------------------------------------------------------------------
module fnvector_test_mod
  use, intrinsic :: iso_c_binding
  use fnvector_complex_mod
  implicit none

contains
  integer(c_int) function check_ans(val, tol, N, sunvec_x) result(failure)

    implicit none
    complex(c_double_complex), value :: val
    real(c_double), value :: tol
    integer(c_long), value :: N
    Type(N_Vector) :: sunvec_x
    Type(FVec), pointer :: x
    integer(c_long) :: i

    x => FN_VGetFVec(sunvec_x)
    failure = 0
    do i = 1,N
       if (abs(x%data(i) - val) > tol)  failure = 1
    end do

  end function check_ans
end module fnvector_test_mod

! ------------------------------------------------------------------
program main

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fnvector_complex_mod
  use fnvector_test_mod

  !======= Declarations =========
  implicit none

  ! local variables
  integer(c_int)  :: fails, i, loc
  integer(c_long), parameter :: N = 1000
  type(N_Vector), pointer :: sU, sV, sW, sX, sY, sZ
  type(FVec), pointer :: U, V, W, X, Y, Z
  complex(c_double_complex) :: Udata(N)
  real(c_double) :: fac
  logical :: failure


  !======= Internals ============

  ! initialize failure total
  fails = 0

  ! create new vectors, using New, Make and Clone routines
  sU => FN_VMake_Complex(N, Udata)
  if (.not. associated(sU)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  U => FN_VGetFVec(sU)

  sV => FN_VNew_Complex(N)
  if (.not. associated(sV)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  V => FN_VGetFVec(sV)

  sW => FN_VNew_Complex(N)
  if (.not. associated(sW)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  W => FN_VGetFVec(sW)

  sX => FN_VNew_Complex(N)
  if (.not. associated(sX)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  X => FN_VGetFVec(sX)

  sY => FN_VNew_Complex(N)
  if (.not. associated(sY)) then
     print *, 'ERROR: sunvec = NULL'
     stop 1
  end if
  Y => FN_VGetFVec(sY)

  call c_f_pointer(FN_VClone_Complex(sU), sZ)
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
  if (FN_VGetLength(sV) /= N) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VGetLength'
     print *, '    ', FN_VGetLength(sV), ' /= ', N
  else
     print *, 'PASSED test -- FN_VGetLength'
  end if

  ! test FN_VConst
  Udata = 0.d0
  call FN_VConst(1.d0, sU)
  if (check_ans(dcmplx(1.d0, 0.d0), 1.d-14, N, sU) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VConst'
  else
     print *, 'PASSED test -- FN_VConst'
  end if

  ! test FN_VLinearSum
  X%data = dcmplx(1.d0, -1.d0)
  Y%data = dcmplx(-2.d0, 2.d0)
  call FN_VLinearSum(1.d0, sX, 1.d0, sY, sY)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sY) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 1a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 1a'
  end if

  X%data = dcmplx(1.d0, -1.d0)
  Y%data = dcmplx(2.d0, -2.d0)
  call FN_VLinearSum(-1.d0, sX, 1.d0, sY, sY)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sY) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 1b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 1b'
  end if

  X%data = dcmplx(2.d0, -2.d0)
  Y%data = dcmplx(-2.d0, 2.d0)
  call FN_VLinearSum(0.5d0, sX, 1.d0, sY, sY)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sY) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 1c'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 1c'
  end if

  X%data = dcmplx(2.d0, -2.d0)
  Y%data = dcmplx(-1.d0, 1.d0)
  call FN_VLinearSum(1.d0, sX, 1.d0, sY, sX)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 2a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 2a'
  end if

  X%data = dcmplx(1.d0, -1.d0)
  Y%data = dcmplx(2.d0, -2.d0)
  call FN_VLinearSum(1.d0, sX, -1.d0, sY, sX)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 2b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 2b'
  end if

  X%data = dcmplx(2.d0, -2.d0)
  Y%data = dcmplx(-0.5d0, 0.5d0)
  call FN_VLinearSum(1.d0, sX, 2.d0, sY, sX)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 2c'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 2c'
  end if

  X%data = dcmplx(-2.d0, 2.d0)
  Y%data = dcmplx(1.d0, -1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(1.d0, sX, 1.d0, sY, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 3'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 3'
  end if

  X%data = dcmplx(2.d0, -2.d0)
  Y%data = dcmplx(1.d0, -1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(1.d0, sX, -1.d0, sY, sZ)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 4a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 4a'
  end if

  X%data = dcmplx(2.d0, -2.d0)
  Y%data = dcmplx(1.d0, -1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(-1.d0, sX, 1.d0, sY, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 4b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 4b'
  end if

  X%data = dcmplx(2.d0, -2.d0)
  Y%data = dcmplx(-0.5d0, 0.5d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(1.d0, sX, 2.d0, sY, sZ)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 5a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 5a'
  end if

  X%data = dcmplx(0.5d0, -0.5d0)
  Y%data = dcmplx(-2.d0, 2.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(2.d0, sX, 1.d0, sY, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 5b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 5b'
  end if

  X%data = dcmplx(-2.d0, 2.d0)
  Y%data = dcmplx(-0.5d0, 0.5d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(-1.d0, sX, 2.d0, sY, sZ)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 6a'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 6a'
  end if

  X%data = dcmplx(0.5d0, -0.5d0)
  Y%data = dcmplx(2.d0, -2.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(2.d0, sX, -1.d0, sY, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 6b'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 6b'
  end if

  X%data = dcmplx(1.d0, -1.d0)
  Y%data = dcmplx(-0.5d0, 0.5d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(2.d0, sX, 2.d0, sY, sZ)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 7'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 7'
  end if

  X%data = dcmplx(0.5d0, -0.5d0)
  Y%data = dcmplx(1.d0, -1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(2.d0, sX, -2.d0, sY, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 8'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 8'
  end if

  X%data = dcmplx(1.d0, -1.d0)
  Y%data = dcmplx(-2.d0, 2.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VLinearSum(2.d0, sX, 0.5d0, sY, sZ)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VLinearSum Case 9'
  else
     print *, 'PASSED test -- FN_VLinearSum Case 9'
  end if

  ! test FN_VProd
  X%data = dcmplx(2.d0, 0.d0)
  Y%data = dcmplx(-0.5d0, 0.0d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VProd(sX, sY, sZ)
  if (check_ans(dcmplx(-1.d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VProd Case 1'
  else
     print *, 'PASSED test -- FN_VProd Case 1'
  end if

  X%data = dcmplx(0.d0, 0.5d0)
  Y%data = dcmplx(-2.0d0, 0.0d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VProd(sX, sY, sZ)
  if (check_ans(dcmplx(0.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VProd Case 2'
  else
     print *, 'PASSED test -- FN_VProd Case 2'
  end if

  X%data = dcmplx(1.d0, 2.d0)
  Y%data = dcmplx(1.0d0, -2.0d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VProd(sX, sY, sZ)
  if (check_ans(dcmplx(5.d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VProd Case 3'
  else
     print *, 'PASSED test -- FN_VProd Case 3'
  end if

  ! test FN_VDiv
  X%data = dcmplx(1.d0, 0.d0)
  Y%data = dcmplx(2.d0, 0.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VDiv(sX, sY, sZ)
  if (check_ans(dcmplx(0.5d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VDiv Case 1'
  else
     print *, 'PASSED test -- FN_VDiv Case 1'
  end if

  X%data = dcmplx(0.d0, 1.d0)
  Y%data = dcmplx(2.d0, 0.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VDiv(sX, sY, sZ)
  if (check_ans(dcmplx(0.d0, 0.5d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VDiv Case 2'
  else
     print *, 'PASSED test -- FN_VDiv Case 2'
  end if

  X%data = dcmplx(4.d0, 2.d0)
  Y%data = dcmplx(1.d0, -1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VDiv(sX, sY, sZ)
  if (check_ans(dcmplx(1.d0, 3.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VDiv Case 3'
  else
     print *, 'PASSED test -- FN_VDiv Case 3'
  end if

  ! test FN_VScale
  X%data = dcmplx(0.5d0, -0.5d0)
  call FN_VScale(2.d0, sX, sX)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sX) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 1'
  else
     print *, 'PASSED test -- FN_VScale Case 1'
  end if

  X%data = dcmplx(-1.d0, 1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VScale(1.d0, sX, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 2'
  else
     print *, 'PASSED test -- FN_VScale Case 2'
  end if

  X%data = dcmplx(-1.d0, 1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VScale(-1.d0, sX, sZ)
  if (check_ans(dcmplx(1.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 3'
  else
     print *, 'PASSED test -- FN_VScale Case 3'
  end if

  X%data = dcmplx(-0.5d0, 0.5d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VScale(2.d0, sX, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VScale Case 4'
  else
     print *, 'PASSED test -- FN_VScale Case 4'
  end if

  ! test FN_VAbs
  X%data = dcmplx(-1.d0, 0.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VAbs(sX, sZ)
  if (check_ans(dcmplx(1.d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VAbs Case 1'
  else
     print *, 'PASSED test -- FN_VAbs Case 1'
  end if

  X%data = dcmplx(1.d0, -0.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VAbs(sX, sZ)
  if (check_ans(dcmplx(1.d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VAbs Case 2'
  else
     print *, 'PASSED test -- FN_VAbs Case 2'
  end if

  X%data = dcmplx(3.d0, -4.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VAbs(sX, sZ)
  if (check_ans(dcmplx(5.d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VAbs Case 3'
  else
     print *, 'PASSED test -- FN_VAbs Case 3'
  end if

  ! test FN_VInv
  X%data = dcmplx(2.d0, 0.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VInv(sX, sZ)
  if (check_ans(dcmplx(0.5d0, 0.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInv Case 1'
  else
     print *, 'PASSED test -- FN_VInv Case 1'
  end if

  X%data = dcmplx(0.d0, 1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VInv(sX, sZ)
  if (check_ans(dcmplx(0.d0, -1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInv Case 2'
  else
     print *, 'PASSED test -- FN_VInv Case 2'
  end if

  ! test FN_VAddConst
  X%data = dcmplx(1.d0, 1.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  call FN_VAddConst(sX, -2.d0, sZ)
  if (check_ans(dcmplx(-1.d0, 1.d0), 1.d-14, N, sZ) /= 0) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VAddConst'
  else
     print *, 'PASSED test -- FN_VAddConst'
  end if

  ! test FN_VMaxNorm
  X%data = dcmplx(-0.5d0, 0.d0)
  X%data(N) = dcmplx(0.d0, -2.d0)
  if (dabs(FN_VMaxNorm(sX) - 2.d0) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VMaxNorm (',FN_VMaxNorm(sX),' /= 2.d0)'
  else
     print *, 'PASSED test -- FN_VMaxNorm'
  end if

  ! test FN_VWrmsNorm
  X%data = dcmplx(-0.5d0, 0.d0)
  Y%data = dcmplx(0.5d0, 0.d0)
  if (dabs(FN_VWrmsNorm(sX,sY) - 0.25d0) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWrmsNorm (',FN_VWrmsNorm(sX,sY),' /= 0.25d0)'
  else
     print *, 'PASSED test -- FN_VWrmsNorm'
  end if

  ! test FN_VWrmsNormMask
  X%data = dcmplx(-0.5d0, 0.d0)
  Y%data = dcmplx(0.5d0, 0.d0)
  Z%data = dcmplx(1.d0, 0.d0)
  Z%data(N) = dcmplx(0.d0, 0.d0)
  fac = dsqrt(1.d0*(N - 1)/N)*0.25d0
  if (dabs(FN_VWrmsNormMask(sX,sY,sZ) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWrmsNormMask (',FN_VWrmsNormMask(sX,sY,sZ),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWrmsNormMask'
  end if

  ! test FN_VMin
  X%data = dcmplx(2.d0, 0.d0)
  X%data(N) = dcmplx(-2.d0, -3.d0)
  if (dabs(FN_VMin(sX) + 2.d0) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VMin (',FN_VMin(sX),' /= -2.d0)'
  else
     print *, 'PASSED test -- FN_VMin'
  end if

  ! test FN_VWL2Norm
  X%data = dcmplx(-0.5d0, 0.d0)
  Y%data = dcmplx(0.5d0, 0.d0)
  fac = dsqrt(1.d0*N)*0.25d0
  if (dabs(FN_VWL2Norm(sX,sY) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWL2Norm (',FN_VWL2Norm(sX,sY),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWL2Norm'
  end if

  ! test FN_VL1Norm
  X%data = dcmplx(0.d0, -1.d0)
  fac = 1.d0*N
  if (dabs(FN_VL1Norm(sX) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VL1Norm (',FN_VL1Norm(sX),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VL1Norm'
  end if

  ! test FN_VInvTest
  X%data = dcmplx(0.5d0, 0.d0)
  Z%data = dcmplx(0.d0, 0.d0)
  failure = (FN_VInvTest(sX, sZ) == 0)
  if ((check_ans(dcmplx(2.d0, 0.d0), 1.d-14, N, sZ) /= 0) .or. failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInvTest Case 1'
  else
     print *, 'PASSED test -- FN_VInvTest Case 1'
  end if

  failure = .false.
  Z%data = dcmplx(0.d0, 0.d0)
  do i = 1,N
     loc = mod(i-1, 2)
     if (loc == 0)  X%data(i) = dcmplx(0.d0, 0.d0)
     if (loc == 1)  X%data(i) = dcmplx(0.5d0, 0.d0)
  end do
  if (FN_VInvTest(sX, sZ) == 1)  failure = .true.
  do i = 1,N
     loc = mod(i-1, 2)
     if ((loc == 0) .and. (Z%data(i) /= dcmplx(0.d0, 0.d0)))  failure = .true.
     if ((loc == 1) .and. (Z%data(i) /= dcmplx(2.d0, 0.d0)))  failure = .true.
  end do
  if (failure) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VInvTest Case 2'
  else
     print *, 'PASSED test -- FN_VInvTest Case 2'
  end if

  ! test FN_VWSqrSumLocal
  X%data = dcmplx(-1.d0, 0.d0)
  Y%data = dcmplx(0.5d0, 0.d0)
  fac = 0.25d0*N
  if (dabs(FN_VWSqrSumLocal(sX,sY) - fac) > 1.d-14) then
     fails = fails + 1
     print *, '>>> FAILED test -- FN_VWSqrSumLocal (',FN_VWSqrSumLocal(sX,sY),' /= ',fac,')'
  else
     print *, 'PASSED test -- FN_VWSqrSumLocal'
  end if


  ! test FN_VWSqrSumMaskLocal
  X%data = dcmplx(-1.d0, 0.d0)
  Y%data = dcmplx(0.5d0, 0.d0)
  Z%data = dcmplx(1.d0, 0.d0)
  Z%data(N) = dcmplx(0.d0, 0.d0)
  fac = 0.25d0*(N-1)
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
