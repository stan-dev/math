! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
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
! This file tests the Fortran 2003 interface to the SUNDIALS
! serial N_Vector implementation.
! -----------------------------------------------------------------

module test_nvector_serial
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_serial_mod
  use test_utilities
  implicit none

  integer(c_long), parameter :: N = 100 ! vector length
  integer(c_int),  parameter :: nv = 3  ! length of vector arrays

contains

  integer function smoke_tests() result(ret)
    implicit none

    integer(c_long)         :: lenrw(1), leniw(1) ! real and int work space size
    integer(c_long)         :: ival               ! integer work value
    type(c_ptr)             :: cptr               ! c_ptr work value
    real(c_double)          :: rval               ! real work value
    real(c_double)          :: xdata(N)           ! vector data array
    real(c_double), pointer :: xptr(:)            ! pointer to vector data array
    real(c_double)          :: nvarr(nv)          ! array of nv constants to go with vector array
    type(N_Vector), pointer :: x, y, z, tmp       ! N_Vectors
    type(c_ptr)             :: xvecs, zvecs       ! C pointer to array of C pointers to N_Vectors

    !===== Setup ====
    x => FN_VMake_Serial(N, xdata, sunctx)
    call FN_VConst(ONE, x)
    y => FN_VClone_Serial(x)
    call FN_VConst(ONE, y)
    z => FN_VClone_Serial(x)
    call FN_VConst(ONE, z)

    xvecs = FN_VCloneVectorArray(nv, x)
    zvecs = FN_VCloneVectorArray(nv, z)
    nvarr = (/ ONE, ONE, ONE /)

    !===== Test =====

    ! test constructors
    tmp => FN_VNewEmpty_Serial(N, sunctx)
    call FN_VDestroy_Serial(tmp)
    tmp => FN_VMake_Serial(N, xdata, sunctx)
    call FN_VDestroy_Serial(tmp)
    tmp => FN_VNew_Serial(N, sunctx)
    call FN_VDestroy_Serial(tmp)
    tmp => FN_VCloneEmpty_Serial(x)
    call FN_VDestroy_Serial(tmp)

    ! test generic vector functions
    ival = FN_VGetVectorID_Serial(x)
    call FN_VSpace_Serial(x, lenrw, leniw)
    xptr => FN_VGetArrayPointer_Serial(x)
    call FN_VSetArrayPointer_Serial(xdata, x)
    cptr = FN_VGetCommunicator(x)
    ival = FN_VGetLength_Serial(x)

    ! test standard vector operations
    call FN_VLinearSum_Serial(ONE, x, ONE, y, z)
    call FN_VConst_Serial(ONE, z)
    call FN_VProd_Serial(x, y, z)
    call FN_VDiv_Serial(x, y, z)
    call FN_VScale_Serial(ONE, x, y)
    call FN_VAbs_Serial(x, y)
    call FN_VInv_Serial(x, z)
    call FN_VAddConst_Serial(x, ONE, z)
    rval = FN_VDotProd_Serial(x, y)
    rval = FN_VMaxNorm_Serial(x)
    rval = FN_VWrmsNorm_Serial(x, y)
    rval = FN_VWrmsNormMask_Serial(x, y, z)
    rval = FN_VMin_Serial(x)
    rval = FN_VWL2Norm_Serial(x, y)
    rval = FN_VL1Norm_Serial(x)
    call FN_VCompare_Serial(ONE, x, y)
    ival = FN_VInvTest_Serial(x, y)
    ival = FN_VConstrMask_Serial(z, x, y)
    rval = FN_VMinQuotient_Serial(x, y)

    ! test fused vector operations
    ival = FN_VLinearCombination_Serial(nv, nvarr, xvecs, x)
    ival = FN_VScaleAddMulti_Serial(nv, nvarr, x, xvecs, zvecs)
    ival = FN_VDotProdMulti_Serial(nv, x, xvecs, nvarr)

    ! test vector array operations
    ival = FN_VLinearSumVectorArray_Serial(nv, ONE, xvecs, ONE, xvecs, zvecs)
    ival = FN_VScaleVectorArray_Serial(nv, nvarr, xvecs, zvecs)
    ival = FN_VConstVectorArray_Serial(nv, ONE, xvecs)
    ival = FN_VWrmsNormVectorArray_Serial(nv, xvecs, xvecs, nvarr)
    ival = FN_VWrmsNormMaskVectorArray_Serial(nv, xvecs, xvecs, x, nvarr)

    !==== Cleanup =====
    call FN_VDestroy_Serial(x)
    call FN_VDestroy_Serial(y)
    call FN_VDestroy_Serial(z)
    call FN_VDestroyVectorArray(xvecs, nv)
    call FN_VDestroyVectorArray(zvecs, nv)

    ret = 0

  end function smoke_tests

  integer function unit_tests() result(fails)
    use test_fnvector
    implicit none

    real(c_double)           :: xdata(N) ! vector data array
    type(N_Vector), pointer  :: x        ! N_Vectors

    !===== Setup ====
    fails = 0

    x => FN_VMake_Serial(N, xdata, sunctx)
    call FN_VConst(ONE, x)

    !==== tests ====
    fails = Test_FN_VMake(x, N, 0)
    fails = Test_FN_VGetArrayPointer(x, N, 0)
    fails = Test_FN_VLinearCombination(x, N, 0)

    !=== cleanup ====
    call FN_VDestroy_Serial(x)

  end function unit_tests

end module


integer(C_INT) function check_ans(ans, X, local_length) result(failure)
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use test_utilities
  implicit none

  real(C_DOUBLE)          :: ans
  type(N_Vector)          :: X
  integer(C_LONG)         :: local_length, i
  real(C_DOUBLE), pointer :: Xdata(:)

  failure = 0

  Xdata => FN_VGetArrayPointer(X)
  do i = 1, local_length
    if (FNEQ(Xdata(i), ans) > 0) then
      failure = failure + 1
    end if
  end do
end function check_ans


logical function has_data(X) result(failure)
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use test_utilities
  implicit none

  type(N_Vector)          :: X
  real(C_DOUBLE), pointer :: xptr(:)

  xptr => FN_VGetArrayPointer(x)
  failure = associated(xptr)
end function has_data


program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_nvector_serial

  !======== Declarations ========
  implicit none
  integer :: fails = 0

  !============== Introduction =============
  print *, 'Serial N_Vector Fortran 2003 interface test'

  call Test_Init(c_null_ptr)

  fails = smoke_tests()
  if (fails /= 0) then
    print *, 'FAILURE: smoke tests failed'
    stop 1
  else
    print *, 'SUCCESS: smoke tests passed'
  end if

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed'
    stop 1
  else
    print *, 'SUCCESS: all unit tests passed'
  end if

  call Test_Finalize()

end program main
