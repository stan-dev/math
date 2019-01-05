! -----------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
! -----------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2014, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! -----------------------------------------------------------------
! This file contains a Fortran module for interfacing directly with
! the SUNDIALS serial NVector using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fnvector_serial_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VNew_Serial
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VNew_Serial(vec_length) &
          bind(C,name='N_VNew_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
     end function FN_VNew_Serial

     ! -----------------------------------------------------------------
     ! N_VNewEmpty_Serial
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VNewEmpty_Serial(vec_length) &
          bind(C,name='N_VNewEmpty_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: vec_length
     end function FN_VNewEmpty_Serial

     ! -----------------------------------------------------------------
     ! N_VMake_Serial
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VMake_Serial(length, v_data) &
          bind(C,name='N_VMake_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: length
       real(c_double)         :: v_data(length)
     end function FN_VMake_Serial

     ! -----------------------------------------------------------------
     ! N_VCloneVectorArray_Serial: NOT INTERFACED
     ! -----------------------------------------------------------------
     
     ! -----------------------------------------------------------------
     ! N_VCloneVectorArrayEmpty_Serial: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Destructors
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VDestroy_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VDestroy_Serial(v) &
          bind(C,name='N_VDestroy_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine FN_VDestroy_Serial

     ! -----------------------------------------------------------------
     ! N_VDestroyVectorArray_Serial: NOT INTERFACED
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Other routines
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VGetLength_Serial
     ! -----------------------------------------------------------------

     integer(c_long) function FN_VGetLength_Serial(v) &
          bind(C,name='N_VGetLength_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end function FN_VGetLength_Serial

     ! -----------------------------------------------------------------
     ! N_VPrint_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VPrint_Serial(v) &
          bind(C,name='N_VPrint_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
     end subroutine FN_VPrint_Serial

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: N_VPrintFile_Serial
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Operations
     ! =================================================================

     ! -----------------------------------------------------------------
     ! N_VGetVectorID_Serial
     ! -----------------------------------------------------------------
     
     integer(c_int) function FN_VGetVectorID_Serial(v) &
          bind(C,name='N_VGetVectorID_Serial')
        use, intrinsic :: iso_c_binding
        implicit none
        type(c_ptr), value :: v
     end function FN_VGetVectorID_Serial
     
     ! -----------------------------------------------------------------
     ! N_VCloneEmpty_Serial
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VCloneEmpty_Serial(w) &
          bind(C,name='N_VCloneEmpty_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function FN_VCloneEmpty_Serial

     ! -----------------------------------------------------------------
     ! N_VClone_Serial
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VClone_Serial(w) &
          bind(C,name='N_VClone_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: w
     end function FN_VClone_Serial

     ! -----------------------------------------------------------------
     ! N_VSpace_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VSpace_Serial(v, lrw, liw) &
          bind(C,name='N_VSpace_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v
       integer(c_long)    :: lrw
       integer(c_long)    :: liw
     end subroutine FN_VSpace_Serial

     ! -----------------------------------------------------------------
     ! N_VGetArrayPointer_Serial
     ! -----------------------------------------------------------------

     type(c_ptr) function FN_VGetArrayPointer_Serial(vec) &
          bind(C,name='N_VGetArrayPointer_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: vec
     end function FN_VGetArrayPointer_Serial

     ! -----------------------------------------------------------------
     ! N_VSetArrayPointer_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VSetArrayPointer_Serial(v_data, v) &
          bind(C,name='N_VSetArrayPointer_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: v_data
       type(c_ptr), value :: v
     end subroutine FN_VSetArrayPointer_Serial

     ! -----------------------------------------------------------------
     ! N_VLinearSum_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VLinearSum_Serial(a, x, b, y, z) &
          bind(C,name='N_VLinearSum_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: a
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: y
       type(c_ptr),    value :: z
     end subroutine FN_VLinearSum_Serial

     ! -----------------------------------------------------------------
     ! N_VConst_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VConst_Serial(c, z) &
          bind(C,name='N_VConst_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: z
     end subroutine FN_VConst_Serial

     ! -----------------------------------------------------------------
     ! N_VProd_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VProd_Serial(x, y, z) &
          bind(C,name='N_VProd_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine FN_VProd_Serial

     ! -----------------------------------------------------------------
     ! N_VDiv_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VDiv_Serial(x, y, z) &
          bind(C,name='N_VDiv_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
       type(c_ptr), value :: z
     end subroutine FN_VDiv_Serial

     ! -----------------------------------------------------------------
     ! N_VScale_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VScale_Serial(c, x, z) &
          bind(C,name='N_VScale_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine FN_VScale_Serial

     ! -----------------------------------------------------------------
     ! N_VAbs_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VAbs_Serial(x, z) &
          bind(C,name='N_VAbs_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine FN_VAbs_Serial

     ! -----------------------------------------------------------------
     ! N_VInv_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VInv_Serial(x, z) &
          bind(C,name='N_VInv_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end subroutine FN_VInv_Serial

     ! -----------------------------------------------------------------
     ! N_VAddConst
     ! -----------------------------------------------------------------

     subroutine FN_VAddConst_Serial(x, b, z) &
          bind(C,name='N_VAddConst_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: x
       real(c_double), value :: b
       type(c_ptr),    value :: z
     end subroutine FN_VAddConst_Serial

     ! -----------------------------------------------------------------
     ! N_VDotProd_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VDotProd_Serial(x, y) &
          bind(C,name='N_VDotProd_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function FN_VDotProd_Serial

     ! -----------------------------------------------------------------
     ! N_VMaxNorm_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMaxNorm_Serial(x) &
          bind(C,name='N_VMaxNorm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VMaxNorm_Serial

     ! -----------------------------------------------------------------
     ! N_VWrmsNorm_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWrmsNorm_Serial(x, w) &
          bind(C,name='N_VWrmsNorm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function FN_VWrmsNorm_Serial

     ! -----------------------------------------------------------------
     ! N_VWrmsNormMask_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWrmsNormMask_Serial(x, w, id) &
          bind(C,name='N_VWrmsNormMask_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
       type(c_ptr), value :: id
     end function FN_VWrmsNormMask_Serial

     ! -----------------------------------------------------------------
     ! N_VMin_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMin_Serial(x) &
          bind(C,name='N_VMin_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VMin_Serial

     ! -----------------------------------------------------------------
     ! N_VWL2Norm_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VWL2Norm_Serial(x, w) &
          bind(C,name='N_VWL2Norm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: w
     end function FN_VWL2Norm_Serial

     ! -----------------------------------------------------------------
     ! N_VL1Norm_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VL1Norm_Serial(x) &
          bind(C,name='N_VL1Norm_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
     end function FN_VL1Norm_Serial

     ! -----------------------------------------------------------------
     ! N_VCompare_Serial
     ! -----------------------------------------------------------------

     subroutine FN_VCompare_Serial(c, x, z) &
          bind(C,name='N_VCompare_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       real(c_double), value :: c
       type(c_ptr),    value :: x
       type(c_ptr),    value :: z
     end subroutine FN_VCompare_Serial

     ! -----------------------------------------------------------------
     ! N_VInvTest_Serial
     ! -----------------------------------------------------------------

     integer(c_int) function FN_VInvTest_Serial(x, z) &
          bind(C,name='N_VInvTest_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: x
       type(c_ptr), value :: z
     end function FN_VInvTest_Serial

     ! -----------------------------------------------------------------
     ! N_VConstrMask_Serial
     ! -----------------------------------------------------------------

     integer(c_int) function FN_VConstrMask_Serial(c, x, m) &
          bind(C,name='N_VConstrMask_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: c
       type(c_ptr), value :: x
       type(c_ptr), value :: m
     end function FN_VConstrMask_Serial

     ! -----------------------------------------------------------------
     ! N_VMinQuotient_Serial
     ! -----------------------------------------------------------------

     real(c_double) function FN_VMinQuotient_Serial(num, denom) &
          bind(C,name='N_VMinQuotient_Serial')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: num
       type(c_ptr), value :: denom
     end function FN_VMinQuotient_Serial

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
  ! FN_VGetData_Serial 
  ! 
  ! Extracts data array from a serial SUNDIALS N_Vector
  ! ----------------------------------------------------------------

  subroutine FN_VGetData_Serial(vec, f_array)

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
    c_array = FN_VGetArrayPointer_Serial(vec)

    ! get vector length
    length = FN_VGetLength_Serial(vec)
    
    ! convert c pointer to f pointer
    call c_f_pointer(c_array, f_array, (/length/))

  end subroutine FN_VGetData_Serial

end module fnvector_serial_mod
