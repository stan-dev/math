! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
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


program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use fnvector_serial_mod

  !======== Declarations ========
  implicit none

  ! constants
  real(c_double) :: ONE = 1.d0
  real(c_double), dimension(100) :: ONES

  ! local variables
  integer(c_int)  :: fails = 0     ! number of test fails
  integer(c_long) :: N = 100       ! vector length
  type(c_ptr)     :: x, y, z       ! NVectors
  integer(c_long) :: lenrw, leniw  ! real and int work space size
  type(c_ptr)     :: cptr          ! a c pointer
  integer(c_long) :: ival          ! integer value
  real(c_double)  :: rval          ! real value
 
  real(c_double), dimension(:), pointer :: rdata  ! fortran real data pointers

  ONES = ONE

  !======= Introduction =========
  print *,'Serial N_Vector Fortran 2003 interface test'

  !===== Calls to interface =====
  
  ! constructors
  x = FN_VNewEmpty_Serial(N)
  if (.not. c_associated(x)) then
    print *,'>>> FAILED - ERROR in FN_VNewEmpty_Serial; halting'
    stop 1
  end if
  call FN_VDestroy_Serial(x)
  
  x = FN_VMake_Serial(N, rdata)
  if (.not. c_associated(x)) then
    print *,'>>> FAILED - ERROR in FN_VMake_Serial; halting'
    stop 1
  end if
  call FN_VDestroy_Serial(x)
  
  x = FN_VNew_Serial(N)
  if (.not. c_associated(x)) then
    print *,'>>> FAILED - ERROR in FN_VNew_Serial; halting'
    stop 1
  end if

  ! zero the vector
  call FN_VConst_Serial(ONE, x)

  ! misc. vector functions
  ival = FN_VGetLength_Serial(x)
  
  ! vector operations 
  ival = FN_VGetVectorID_Serial(x)

  y = FN_VCloneEmpty_Serial(x)
  if (.not. c_associated(y)) then
    print *,'>>> FAILED - ERROR in FN_VCloneEmpty_Serial; halting'
    stop 1
  end if
  call FN_VDestroy_Serial(y)

  y = FN_VClone_Serial(x)
  if (.not. c_associated(y)) then
    print *,'>>> FAILED - ERROR in FN_VClone_Serial; halting'
    stop 1
  end if
  
  z = FN_VClone_Serial(x)
  if (.not. c_associated(z)) then
    print *,'>>> FAILED - ERROR in FN_VClone_Serial; halting'
    stop 1
  end if

  call FN_VSpace_Serial(x, lenrw, leniw)
  cptr = FN_VGetArrayPointer_Serial(x)
  call FN_VSetArrayPointer_Serial(cptr, x)
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

  !======= Cleanup ===========
  call FN_VDestroy_Serial(x)
  call FN_VDestroy_Serial(y)
  call FN_VDestroy_Serial(z)

  if (fails == 0) then
    print *,'    SUCCESS - all tests passed'
  else
    print *,'    FAILURE - ', fails, ' tests failed'
    stop 1
  end if

end program main

