! -----------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
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
! the SUNDIALS band matrix using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunmatrix_band_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNBandMatrix(N, mu, ml, smu) &
          bind(C,name='SUNBandMatrix')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: N
       integer(c_long), value :: mu
       integer(c_long), value :: ml
       integer(c_long), value :: smu
     end function FSUNBandMatrix

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNMatDestroy_Band(A) &
          bind(C,name='SUNMatDestroy_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end subroutine FSUNMatDestroy_Band

     ! =================================================================
     ! Other routines
     ! =================================================================

     ! -----------------------------------------------------------------
     ! NOT INTERFACED SUNBandMatrix_Print
     ! -----------------------------------------------------------------

     integer(c_long) function FSUNBandMatrix_Rows(A) &
          bind(C,name='SUNBandMatrix_Rows')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_Rows

     integer(c_long) function FSUNBandMatrix_Columns(A) &
          bind(C,name='SUNBandMatrix_Columns')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_Columns

     integer(c_long) function FSUNBandMatrix_LowerBandwidth(A) &
          bind(C,name='SUNBandMatrix_LowerBandwidth')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_LowerBandwidth

     integer(c_long) function FSUNBandMatrix_UpperBandwidth(A) &
          bind(C,name='SUNBandMatrix_UpperBandwidth')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_UpperBandwidth

     integer(c_long) function FSUNBandMatrix_StoredUpperBandwidth(A) &
          bind(C,name='SUNBandMatrix_StoredUpperBandwidth')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_StoredUpperBandwidth

     integer(c_long) function FSUNBandMatrix_LDim(A) &
          bind(C,name='SUNBandMatrix_LDim')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_LDim

     type(c_ptr) function FSUNBandMatrix_Data(A) &
          bind(C,name='SUNBandMatrix_Data')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNBandMatrix_Data

     ! -----------------------------------------------------------------
     ! NOT INTERFACED SUNBandMatrix_Cols
     ! -----------------------------------------------------------------

     type(c_ptr) function FSUNBandMatrix_Column(A, j) &
          bind(C,name='SUNBandMatrix_Column')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value     :: A
       integer(c_long), value :: j
     end function FSUNBandMatrix_Column

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNMatGetID_Band(A) &
          bind(C,name='SUNMatGetID_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNMatGetID_Band

     type(c_ptr) function FSUNMatClone_Band(A) &
          bind(C,name='SUNMatClone_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNMatClone_Band

     integer(c_int) function FSUNMatZero_Band(A) &
          bind(C,name='SUNMatZero_Band')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
     end function FSUNMatZero_Band

     integer(c_int) function FSUNMatCopy_Band(A, B) &
          bind(C,name='SUNMatCopy_Band')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function FSUNMatCopy_Band

     integer(c_int) function FSUNMatScaleAdd_Band(c, A, B) &
          bind(C,name='SUNMatScaleAdd_Band')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: c
       type(c_ptr),    value :: A
       type(c_ptr),    value :: B
     end function FSUNMatScaleAdd_Band

     integer(c_int) function FSUNMatScaleAddI_Band(c, A) &
          bind(C,name='SUNMatScaleAddI_Band')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: c
       type(c_ptr),    value :: A
     end function FSUNMatScaleAddI_Band

     integer(c_int) function FSUNMatMatvec_Band(A, x, y) &
          bind(C,name='SUNMatMatvec_Band')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function FSUNMatMatvec_Band

     integer(c_int) function FSUNMatSpace_Band(A, lenrw, leniw) &
          bind(C,name='SUNMatSpace_Band')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       integer(c_long)    :: lenrw
       integer(c_long)    :: leniw
     end function FSUNMatSpace_Band

  end interface

contains

  ! ================================================================
  ! Helpful routines
  ! ================================================================

  ! ----------------------------------------------------------------
  ! FSUNMatGetData_Band
  !
  ! Extracts data array from a SUNDIALS Band Matrix
  ! ----------------------------------------------------------------

  subroutine FSUNMatGetData_Band(A, f_array)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr)             :: A
    real(c_double), pointer :: f_array(:)

    ! internal variables
    type(c_ptr)     :: c_array
    integer(c_long) :: ldim

    !======= Internals ============
    
    ! get data pointer from N_Vector
    c_array = FSUNBandMatrix_Data(A)

    ! get length of the data array
    ldim = FSUNBandMatrix_LDim(A)

    ! convert 1D data array
    call c_f_pointer(c_array, f_array, (/ldim/))

  end subroutine FSUNMatGetData_Band

end module fsunmatrix_band_mod
