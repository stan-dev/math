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
! This file contains a Fortran module for interfacing directly with
! the SUNDIALS Pthreads NVector using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fnvector_pthreads_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VNew_Pthreads
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VNew_Pthreads(vec_length, num_threads) &
          bind(C,name='N_VNew_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
       integer(c_int),  value :: num_threads
     end function FN_VNew_Pthreads

     ! -----------------------------------------------------------------
     ! N_VNewEmpty_Pthreads
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VNewEmpty_Pthreads(vec_length, num_threads) &
          bind(C,name='N_VNewEmpty_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
       integer(c_int),  value :: num_threads
     end function FN_VNewEmpty_Pthreads

     ! -----------------------------------------------------------------
     ! N_VMake_Pthreads
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VMake_Pthreads(length, v_data, num_threads) &
          bind(C,name='N_VMake_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: length
       real(c_double)         :: v_data(length)
       integer(c_int),  value :: num_threads
     end function FN_VMake_Pthreads

     ! -----------------------------------------------------------------
     ! N_VCloneVectorArray_Pthreads: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! -----------------------------------------------------------------
     ! N_VCloneVectorArrayEmpty_Pthreads: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Destructors
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VDestroy_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VDestroy_Pthreads(v) &
          bind(C,name='N_VDestroy_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine FN_VDestroy_Pthreads

     ! -----------------------------------------------------------------
     ! N_VDestroyVectorArray_Pthreads: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Other routines
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VGetLength_Pthreads
     ! -----------------------------------------------------------------

     integer(c_long) function FN_VGetLength_Pthreads(v) &
          bind(C,name='N_VGetLength_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end function FN_VGetLength_Pthreads

     ! -----------------------------------------------------------------
     ! N_VPrint_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VPrint_Pthreads(v) &
          bind(C,name='N_VPrint_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine FN_VPrint_Pthreads

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: N_VPrintFile_Pthreads
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Operations
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VGetVectorID_Pthreads
     ! -----------------------------------------------------------------
     
     integer(c_int) function FN_VGetVectorID_Pthreads(v) &
          bind(C,name='N_VGetVectorID_Pthreads')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), value :: v
     end function FN_VGetVectorID_Pthreads

     ! -----------------------------------------------------------------
     ! N_VCloneEmpty_Pthreads
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VCloneEmpty_Pthreads(w) &
          bind(C,name='N_VCloneEmpty_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function FN_VCloneEmpty_Pthreads

     ! -----------------------------------------------------------------
     ! N_VClone_Pthreads
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VClone_Pthreads(w) &
          bind(C,name='N_VClone_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function FN_VClone_Pthreads

     ! -----------------------------------------------------------------
     ! N_VSpace_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VSpace_Pthreads(v, lrw, liw) &
          bind(C,name='N_VSpace_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
       integer(c_long)    :: lrw
       integer(c_long)    :: liw
     end subroutine FN_VSpace_Pthreads

     ! -----------------------------------------------------------------
     ! N_VGetArrayPointer_Pthreads
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VGetArrayPointer_Pthreads(vec) &
          bind(C,name='N_VGetArrayPointer_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vec
     end function FN_VGetArrayPointer_Pthreads

     ! -----------------------------------------------------------------
     ! N_VSetArrayPointer_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VSetArrayPointer_Pthreads(v_data, v) &
          bind(C,name='N_VSetArrayPointer_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v_data
       type(c_ptr), value :: v
     end subroutine FN_VSetArrayPointer_Pthreads

     ! -----------------------------------------------------------------
     ! N_VLinearSum_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VLinearSum_Pthreads(a, x, b, y, z) &
          bind(C,name='N_VLinearSum_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: a
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: y
       type(c_ptr),    value :: z
     end subroutine FN_VLinearSum_Pthreads

     ! -----------------------------------------------------------------
     ! N_VConst_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VConst_Pthreads(c, z) &
          bind(C,name='N_VConst_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: z
     end subroutine FN_VConst_Pthreads

     ! -----------------------------------------------------------------
     ! N_VProd_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VProd_Pthreads(x, y, z) &
          bind(C,name='N_VProd_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine FN_VProd_Pthreads

     ! -----------------------------------------------------------------
     ! N_VDiv_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VDiv_Pthreads(x, y, z) &
          bind(C,name='N_VDiv_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine FN_VDiv_Pthreads

     ! -----------------------------------------------------------------
     ! N_VScale_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VScale_Pthreads(c, x, z) &
          bind(C,name='N_VScale_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine FN_VScale_Pthreads

     ! -----------------------------------------------------------------
     ! N_VAbs_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VAbs_Pthreads(x, z) &
          bind(C,name='N_VAbs_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine FN_VAbs_Pthreads

     ! -----------------------------------------------------------------
     ! N_VInv_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VInv_Pthreads(x, z) &
          bind(C,name='N_VInv_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine FN_VInv_Pthreads

     ! -----------------------------------------------------------------
     ! N_VAddConst
     ! -----------------------------------------------------------------

     subroutine FN_VAddConst_Pthreads(x, b, z) &
          bind(C,name='N_VAddConst_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: z
     end subroutine FN_VAddConst_Pthreads

     ! -----------------------------------------------------------------
     ! N_VDotProd_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VDotProd_Pthreads(x, y) &
          bind(C,name='N_VDotProd_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function FN_VDotProd_Pthreads

     ! -----------------------------------------------------------------
     ! N_VMaxNorm_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMaxNorm_Pthreads(x) &
          bind(C,name='N_VMaxNorm_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VMaxNorm_Pthreads

     ! -----------------------------------------------------------------
     ! N_VWrmsNorm_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWrmsNorm_Pthreads(x, w) &
          bind(C,name='N_VWrmsNorm_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function FN_VWrmsNorm_Pthreads

     ! -----------------------------------------------------------------
     ! N_VWrmsNormMask_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWrmsNormMask_Pthreads(x, w, id) &
          bind(C,name='N_VWrmsNormMask_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
       type(c_ptr), value :: id
     end function FN_VWrmsNormMask_Pthreads

     ! -----------------------------------------------------------------
     ! N_VMin_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMin_Pthreads(x) &
          bind(C,name='N_VMin_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VMin_Pthreads

     ! -----------------------------------------------------------------
     ! N_VWL2Norm_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWL2Norm_Pthreads(x, w) &
          bind(C,name='N_VWL2Norm_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function FN_VWL2Norm_Pthreads

     ! -----------------------------------------------------------------
     ! N_VL1Norm_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VL1Norm_Pthreads(x) &
          bind(C,name='N_VL1Norm_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VL1Norm_Pthreads

     ! -----------------------------------------------------------------
     ! N_VCompare_Pthreads
     ! -----------------------------------------------------------------

     subroutine FN_VCompare_Pthreads(c, x, z) &
          bind(C,name='N_VCompare_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine FN_VCompare_Pthreads

     ! -----------------------------------------------------------------
     ! N_VInvTest_Pthreads
     ! -----------------------------------------------------------------

     integer(c_int) function FN_VInvTest_Pthreads(x, z) &
          bind(C,name='N_VInvTest_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end function FN_VInvTest_Pthreads

     ! -----------------------------------------------------------------
     ! N_VConstrMask_Pthreads
     ! -----------------------------------------------------------------

     integer(c_int) function FN_VConstrMask_Pthreads(c, x, m) &
          bind(C,name='N_VConstrMask_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: c
       type(c_ptr), value :: x
       type(c_ptr), value :: m
     end function FN_VConstrMask_Pthreads

     ! -----------------------------------------------------------------
     ! N_VMinQuotient_Pthreads
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMinQuotient_Pthreads(num, denom) &
          bind(C,name='N_VMinQuotient_Pthreads')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: num
       type(c_ptr), value :: denom
     end function FN_VMinQuotient_Pthreads

     ! ================================================================
     ! Fused vector operations: NOT INTERFACED
     ! ================================================================

     ! ================================================================
     ! Vector array operations: NOT INTERFACED
     ! ================================================================
     
   end interface

contains

  ! ================================================================
  ! Helpful routines
  ! ================================================================

  ! ----------------------------------------------------------------
  ! FN_VGetData_Pthreads 
  ! 
  ! Extracts data array from a Pthreads SUNDIALS N_Vector
  ! ----------------------------------------------------------------

  subroutine FN_VGetData_Pthreads(vec, f_array)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr)             :: vec
    integer(c_long)         :: length
    real(c_double), pointer :: f_array(:)

    ! C pointer for N_Vector interal data array
    type(c_ptr) :: c_array

    !======= Internals ============

    ! get data pointer from N_Vector
    c_array = FN_VGetArrayPointer_Pthreads(vec)

    ! get vector length
    length = FN_VGetLength_Pthreads(vec)
    
    ! convert c pointer to f pointer
    call c_f_pointer(c_array, f_array, (/length/))

  end subroutine FN_VGetData_Pthreads

end module fnvector_pthreads_mod
