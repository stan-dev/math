! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! Acknowledgements: These testing routines are based on
! test_nvector.c written by David Gardner and Slaven Peles @ LLNL.
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2022, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------
! These test functions are designed to check the fortran interface
! to an NVECTOR module implementation. It does not test every
! function. It tests the N_VMake constructor, one standard vector
! operation (N_VConst), N_VGetArrayPointer, and one fused operation.
! -----------------------------------------------------------------

module test_fnvector
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fsundials_types_mod
  use test_utilities
  implicit none

  integer(C_INT), external :: check_ans
  logical, external        :: has_data

contains


integer(C_INT) function Test_FN_VMake(X, local_length, myid) &
    result(failure)
  implicit none

  type(N_Vector)  :: X
  integer(C_LONG) :: local_length
  integer(C_INT)  :: myid

  if (.not. has_data(X)) then
    print *, '(I4)', '>>> FAILED test -- FN_VMake, Proc ', myid
    print *, '    vector data is not associated'
    failure = 1
    return
  end if

  if (myid == 0) then
    print *, 'PASSED test -- FN_VMake'
  end if

  failure = 0
end function Test_FN_VMake


!! ----------------------------------------------------------------------
!! NOTE: This routine depends on FN_VConst to check vector data.
!! ----------------------------------------------------------------------
integer(C_INT) function Test_FN_VGetArrayPointer(W, local_length, myid) &
    result(failure)
  implicit none

  type(N_Vector)  :: W
  integer(C_LONG) :: local_length
  integer(C_INT)  :: myid

  ! check vector data
  if (.not. has_data(W)) then
    print *, '>>> FAILED test -- FN_VGetArrayPointer, Proc ', myid
    print *, '    Vector data == NULL \n\n'
    failure = 1
    return;
  end if

  call FN_VConst(NEG_HALF, W)
  failure = check_ans(NEG_HALF, W, local_length)

  if (failure > 0) then
    print *, '(I2)', '>>> FAILED test -- FN_VGetArrayPointer, Proc ', myid
    print *, '    Failed FN_VConst check \n\n'
    failure = 1
    return
  end if

  if (myid == 0) then
    print *, 'PASSED test -- FN_VConst'
    print *, 'PASSED test -- FN_VGetArrayPointer'
  end if

  failure = 0
end function Test_FN_VGetArrayPointer


integer(C_INT) function Test_FN_VLinearCombination(X, local_length, myid) &
    result(failure)

  type(N_Vector)          :: X
  integer(C_LONG)         :: local_length
  integer(C_INT)          :: myid, ierr
  type(N_Vector), pointer :: Y1, Y2, Y3
  type(c_ptr), target     :: V(3)
  type(c_ptr)             :: Vptr
  real(C_DOUBLE)          :: c(3)

  failure = 0

  ! create vectors for testing
  Y1 => FN_VClone(X)
  Y2 => FN_VClone(X)
  Y3 => FN_VClone(X)

  ! set vectors in vector array
  V = (/c_loc(Y1), c_loc(Y2), c_loc(Y3)/)
  Vptr = c_loc(V)

  ! initialize c values
  c = ZERO

  !
  ! Case 1a: V[0] = a V[0], FN_VScale
  !

  ! fill vector data
  call FN_VConst(TWO, Y1)

  ! set scaling factors
  c =  HALF

  ierr = FN_VLinearCombination(1, c, Vptr, Y1)

  ! Y1 should be vector of +1
  if (ierr == 0) then
    failure = check_ans(ONE, Y1, local_length)
  else
    failure = 1
  end if

  if (failure > 0) then
    print *, '(I4)', '>>> FAILED test -- FN_VLinearCombination Case 1a, Proc ', myid
  else if (myid == 0) then
    print *, 'PASSED test -- FN_VLinearCombination Case 1a'
  end if

  !
  ! Case 3a: V[0] = V[0] + b V[1] + c V[2]
  !

  call FN_VConst(TWO, Y1)
  call FN_VConst(NEG_TWO, Y2)
  call FN_VConst(NEG_ONE, Y3)

  c(1) = ONE
  c(2) = HALF
  c(3) = NEG_TWO

  ierr = FN_VLinearCombination(3, c, Vptr, Y1)

  ! Y1 should be vector of +3
  if (ierr == 0) then
    failure = check_ans(TWO+ONE, Y1, local_length)
  else
    failure = 1
  end if

  if (failure > 0) then
    print *, '(I4)', '>>> FAILED test -- FN_VLinearCombination Case 3a, Proc ', myid
  else if (myid == 0) then
    print *, 'PASSED test -- FN_VLinearCombination Case 3a'
  end if

  ! Free vectors
  call FN_VDestroy(Y1);
  call FN_VDestroy(Y2);
  call FN_VDestroy(Y3);

end function Test_FN_VLinearCombination

end module
