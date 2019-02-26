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
! the SUNDIALS OpenMP NVector using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fnvector_openmp_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VNew_OpenMP
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VNew_OpenMP(vec_length, num_threads) &
          bind(C,name='N_VNew_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
       integer(c_int),  value :: num_threads
     end function FN_VNew_OpenMP

     ! -----------------------------------------------------------------
     ! N_VNewEmpty_OpenMP
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VNewEmpty_OpenMP(vec_length, num_threads) &
          bind(C,name='N_VNewEmpty_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
       integer(c_int),  value :: num_threads
     end function FN_VNewEmpty_OpenMP

     ! -----------------------------------------------------------------
     ! N_VMake_OpenMP
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VMake_OpenMP(length, v_data, num_threads) &
          bind(C,name='N_VMake_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: length
       real(c_double)         :: v_data(length)
       integer(c_int),  value :: num_threads
     end function FN_VMake_OpenMP

     ! -----------------------------------------------------------------
     ! N_VCloneVectorArray_OpenMP: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! -----------------------------------------------------------------
     ! N_VCloneVectorArrayEmpty_OpenMP: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Destructors
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VDestroy_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VDestroy_OpenMP(v) &
          bind(C,name='N_VDestroy_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine FN_VDestroy_OpenMP

     ! -----------------------------------------------------------------
     ! N_VDestroyVectorArray_OpenMP: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Other routines
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VGetLength_OpenMP
     ! -----------------------------------------------------------------

     integer(c_long) function FN_VGetLength_OpenMP(v) &
          bind(C,name='N_VGetLength_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end function FN_VGetLength_OpenMP

     ! -----------------------------------------------------------------
     ! N_VPrint_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VPrint_OpenMP(v) &
          bind(C,name='N_VPrint_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine FN_VPrint_OpenMP

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: N_VPrintFile_OpenMP
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Operations
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VGetVectorID_OpenMP
     ! -----------------------------------------------------------------
     
     integer(c_int) function FN_VGetVectorID_OpenMP(v) &
          bind(C,name='N_VGetVectorID_OpenMP')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), value :: v
     end function FN_VGetVectorID_OpenMP

     ! -----------------------------------------------------------------
     ! N_VCloneEmpty_OpenMP
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VCloneEmpty_OpenMP(w) &
          bind(C,name='N_VCloneEmpty_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function FN_VCloneEmpty_OpenMP

     ! -----------------------------------------------------------------
     ! N_VClone_OpenMP
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VClone_OpenMP(w) &
          bind(C,name='N_VClone_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function FN_VClone_OpenMP

     ! -----------------------------------------------------------------
     ! N_VSpace_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VSpace_OpenMP(v, lrw, liw) &
          bind(C,name='N_VSpace_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
       integer(c_long)    :: lrw
       integer(c_long)    :: liw
     end subroutine FN_VSpace_OpenMP

     ! -----------------------------------------------------------------
     ! N_VGetArrayPointer_OpenMP
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VGetArrayPointer_OpenMP(vec) &
          bind(C,name='N_VGetArrayPointer_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vec
     end function FN_VGetArrayPointer_OpenMP

     ! -----------------------------------------------------------------
     ! N_VSetArrayPointer_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VSetArrayPointer_OpenMP(v_data, v) &
          bind(C,name='N_VSetArrayPointer_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v_data
       type(c_ptr), value :: v
     end subroutine FN_VSetArrayPointer_OpenMP

     ! -----------------------------------------------------------------
     ! N_VLinearSum_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VLinearSum_OpenMP(a, x, b, y, z) &
          bind(C,name='N_VLinearSum_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: a
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: y
       type(c_ptr),    value :: z
     end subroutine FN_VLinearSum_OpenMP

     ! -----------------------------------------------------------------
     ! N_VConst_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VConst_OpenMP(c, z) &
          bind(C,name='N_VConst_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: z
     end subroutine FN_VConst_OpenMP

     ! -----------------------------------------------------------------
     ! N_VProd_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VProd_OpenMP(x, y, z) &
          bind(C,name='N_VProd_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine FN_VProd_OpenMP

     ! -----------------------------------------------------------------
     ! N_VDiv_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VDiv_OpenMP(x, y, z) &
          bind(C,name='N_VDiv_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine FN_VDiv_OpenMP

     ! -----------------------------------------------------------------
     ! N_VScale_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VScale_OpenMP(c, x, z) &
          bind(C,name='N_VScale_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine FN_VScale_OpenMP

     ! -----------------------------------------------------------------
     ! N_VAbs_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VAbs_OpenMP(x, z) &
          bind(C,name='N_VAbs_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine FN_VAbs_OpenMP

     ! -----------------------------------------------------------------
     ! N_VInv_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VInv_OpenMP(x, z) &
          bind(C,name='N_VInv_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine FN_VInv_OpenMP

     ! -----------------------------------------------------------------
     ! N_VAddConst
     ! -----------------------------------------------------------------

     subroutine FN_VAddConst_OpenMP(x, b, z) &
          bind(C,name='N_VAddConst_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: z
     end subroutine FN_VAddConst_OpenMP

     ! -----------------------------------------------------------------
     ! N_VDotProd_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VDotProd_OpenMP(x, y) &
          bind(C,name='N_VDotProd_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function FN_VDotProd_OpenMP

     ! -----------------------------------------------------------------
     ! N_VMaxNorm_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMaxNorm_OpenMP(x) &
          bind(C,name='N_VMaxNorm_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VMaxNorm_OpenMP

     ! -----------------------------------------------------------------
     ! N_VWrmsNorm_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWrmsNorm_OpenMP(x, w) &
          bind(C,name='N_VWrmsNorm_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function FN_VWrmsNorm_OpenMP

     ! -----------------------------------------------------------------
     ! N_VWrmsNormMask_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWrmsNormMask_OpenMP(x, w, id) &
          bind(C,name='N_VWrmsNormMask_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
       type(c_ptr), value :: id
     end function FN_VWrmsNormMask_OpenMP

     ! -----------------------------------------------------------------
     ! N_VMin_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMin_OpenMP(x) &
          bind(C,name='N_VMin_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VMin_OpenMP

     ! -----------------------------------------------------------------
     ! N_VWL2Norm_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWL2Norm_OpenMP(x, w) &
          bind(C,name='N_VWL2Norm_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function FN_VWL2Norm_OpenMP

     ! -----------------------------------------------------------------
     ! N_VL1Norm_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VL1Norm_OpenMP(x) &
          bind(C,name='N_VL1Norm_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VL1Norm_OpenMP

     ! -----------------------------------------------------------------
     ! N_VCompare_OpenMP
     ! -----------------------------------------------------------------

     subroutine FN_VCompare_OpenMP(c, x, z) &
          bind(C,name='N_VCompare_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine FN_VCompare_OpenMP

     ! -----------------------------------------------------------------
     ! N_VInvTest_OpenMP
     ! -----------------------------------------------------------------

     integer(c_int) function FN_VInvTest_OpenMP(x, z) &
          bind(C,name='N_VInvTest_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end function FN_VInvTest_OpenMP

     ! -----------------------------------------------------------------
     ! N_VConstrMask_OpenMP
     ! -----------------------------------------------------------------

     integer(c_int) function FN_VConstrMask_OpenMP(c, x, m) &
          bind(C,name='N_VConstrMask_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: c
       type(c_ptr), value :: x
       type(c_ptr), value :: m
     end function FN_VConstrMask_OpenMP

     ! -----------------------------------------------------------------
     ! N_VMinQuotient_OpenMP
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMinQuotient_OpenMP(num, denom) &
          bind(C,name='N_VMinQuotient_OpenMP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: num
       type(c_ptr), value :: denom
     end function FN_VMinQuotient_OpenMP

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
  ! FN_VGetData_OpenMP 
  ! 
  ! Extracts data array from a OpenMP SUNDIALS N_Vector
  ! ----------------------------------------------------------------

  subroutine FN_VGetData_OpenMP(vec, f_array)

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
    c_array = FN_VGetArrayPointer_OpenMP(vec)

    ! get vector length
    length = FN_VGetLength_OpenMP(vec)
    
    ! convert c pointer to f pointer
    call c_f_pointer(c_array, f_array, (/length/))

  end subroutine FN_VGetData_OpenMP

end module fnvector_openmp_mod
