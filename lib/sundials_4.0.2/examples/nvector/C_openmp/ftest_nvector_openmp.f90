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
! OpenMP N_Vector implementation.
! -----------------------------------------------------------------


program main

  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use fnvector_openmp_mod

  !======== Declarations ========
  implicit none

  ! constants
  real(c_double) :: ONE = 1.d0
  real(c_double), dimension(100) :: ONES

  ! local variables
  integer(c_int)  :: fails = 0     ! number of test fails
  integer(c_long) :: N = 100       ! vector length
  integer(c_int)  :: nthreads = 1  ! number of threads
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
  x = FN_VNewEmpty_OpenMP(N, nthreads)
  if (.not. c_associated(x)) then
    print *,'>>> FAILED - ERROR in FN_VNewEmpty_OpenMP; halting'
    stop 1
  end if
  call FN_VDestroy_OpenMP(x)
  
  x = FN_VMake_OpenMP(N, rdata, nthreads)
  if (.not. c_associated(x)) then
    print *,'>>> FAILED - ERROR in FN_VMake_OpenMP; halting'
    stop 1
  end if
  call FN_VDestroy_OpenMP(x)
  
  x = FN_VNew_OpenMP(N, nthreads)
  if (.not. c_associated(x)) then
    print *,'>>> FAILED - ERROR in FN_VNew_OpenMP; halting'
    stop 1
  end if

  ! zero the vector
  call FN_VConst_OpenMP(ONE, x)

  ! misc. vector functions
  ival = FN_VGetLength_OpenMP(x)
  
  ! vector operations 
  ival = FN_VGetVectorID_OpenMP(x)
  
  y = FN_VCloneEmpty_OpenMP(x)
  if (.not. c_associated(y)) then
    print *,'>>> FAILED - ERROR in FN_VCloneEmpty_OpenMP; halting'
    stop 1
  end if
  call FN_VDestroy_OpenMP(y)

  y = FN_VClone_OpenMP(x)
  if (.not. c_associated(y)) then
    print *,'>>> FAILED - ERROR in FN_VClone_OpenMP; halting'
    stop 1
  end if
  
  z = FN_VClone_OpenMP(x)
  if (.not. c_associated(z)) then
    print *,'>>> FAILED - ERROR in FN_VClone_OpenMP; halting'
    stop 1
  end if

  call FN_VSpace_OpenMP(x, lenrw, leniw)
  cptr = FN_VGetArrayPointer_OpenMP(x)
  call FN_VSetArrayPointer_OpenMP(cptr, x)
  call FN_VLinearSum_OpenMP(ONE, x, ONE, y, z)
  call FN_VConst_OpenMP(ONE, z)
  call FN_VProd_OpenMP(x, y, z)
  call FN_VDiv_OpenMP(x, y, z)
  call FN_VScale_OpenMP(ONE, x, y)
  call FN_VAbs_OpenMP(x, y)
  call FN_VInv_OpenMP(x, z)
  call FN_VAddConst_OpenMP(x, ONE, z)
  rval = FN_VDotProd_OpenMP(x, y)
  rval = FN_VMaxNorm_OpenMP(x)
  rval = FN_VWrmsNorm_OpenMP(x, y)
  rval = FN_VWrmsNormMask_OpenMP(x, y, z)
  rval = FN_VMin_OpenMP(x)
  rval = FN_VWL2Norm_OpenMP(x, y)
  rval = FN_VL1Norm_OpenMP(x)
  call FN_VCompare_OpenMP(ONE, x, y)
  ival = FN_VInvTest_OpenMP(x, y)
  ival = FN_VConstrMask_OpenMP(z, x, y)
  rval = FN_VMinQuotient_OpenMP(x, y)

  !======= Cleanup ===========
  call FN_VDestroy_OpenMP(x)
  call FN_VDestroy_OpenMP(y)
  call FN_VDestroy_OpenMP(z)

  if (fails == 0) then
    print *,'    SUCCESS - all tests passed'
  else
    print *,'    FAILURE - ', fails, ' tests failed'
    stop 1
  end if

end program main

