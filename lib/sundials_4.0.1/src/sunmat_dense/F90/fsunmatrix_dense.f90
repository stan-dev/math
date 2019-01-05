! -----------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
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
! the SUNDIALS dense matrix using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunmatrix_dense_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNDenseMatrix(M, N) &
          bind(C,name='SUNDenseMatrix')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: M
       integer(c_long), value :: N
     end function FSUNDenseMatrix

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNMatDestroy_Dense(A) &
          bind(C,name='SUNMatDestroy_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end subroutine FSUNMatDestroy_Dense

     ! =================================================================
     ! Other routines
     ! =================================================================

     ! -----------------------------------------------------------------
     ! NOT INTERFACED SUNDenseMatrix_Print
     ! -----------------------------------------------------------------

     integer(c_long) function FSUNDenseMatrix_Rows(A) &
          bind(C,name='SUNDenseMatrix_Rows')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNDenseMatrix_Rows

     integer(c_long) function FSUNDenseMatrix_Columns(A) &
          bind(C,name='SUNDenseMatrix_Columns')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNDenseMatrix_Columns

     integer(c_long) function FSUNDenseMatrix_LData(A) &
          bind(C,name='SUNDenseMatrix_LData')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNDenseMatrix_LData

     type(c_ptr) function FSUNDenseMatrix_Data(A) &
          bind(C,name='SUNDenseMatrix_Data')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNDenseMatrix_Data

     ! -----------------------------------------------------------------
     ! NOT INTERFACED SUNDenseMatrix_Cols
     ! -----------------------------------------------------------------

     type(c_ptr) function FSUNDenseMatrix_Column(A, j) &
          bind(C,name='SUNDenseMatrix_Column')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value     :: A
       integer(c_long), value :: j
     end function FSUNDenseMatrix_Column

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNMatGetID_Dense(A) &
          bind(C,name='SUNMatGetID_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNMatGetID_Dense

     type(c_ptr) function FSUNMatClone_Dense(A) &
          bind(C,name='SUNMatClone_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNMatClone_Dense

     integer(c_int) function FSUNMatZero_Dense(A) &
          bind(C,name='SUNMatZero_Dense')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
     end function FSUNMatZero_Dense

     integer(c_int) function FSUNMatCopy_Dense(A, B) &
          bind(C,name='SUNMatCopy_Dense')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function FSUNMatCopy_Dense

     integer(c_int) function FSUNMatScaleAdd_Dense(c, A, B) &
          bind(C,name='SUNMatScaleAdd_Dense')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: c
       type(c_ptr),    value :: A
       type(c_ptr),    value :: B
     end function FSUNMatScaleAdd_Dense

     integer(c_int) function FSUNMatScaleAddI_Dense(c, A) &
          bind(C,name='SUNMatScaleAddI_Dense')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: c
       type(c_ptr),    value :: A
     end function FSUNMatScaleAddI_Dense

     integer(c_int) function FSUNMatMatvec_Dense(A, x, y) &
          bind(C,name='SUNMatMatvec_Dense')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function FSUNMatMatvec_Dense

     integer(c_int) function FSUNMatSpace_Dense(A, lenrw, leniw) &
          bind(C,name='SUNMatSpace_Dense')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       integer(c_long)    :: lenrw
       integer(c_long)    :: leniw
     end function FSUNMatSpace_Dense

  end interface

contains

  ! ================================================================
  ! Helpful routines
  ! ================================================================

  ! ----------------------------------------------------------------
  ! FSUNMatGetData_Dense
  !
  ! Extracts data array from a SUNDIALS Dense Matrix
  ! ----------------------------------------------------------------

  subroutine FSUNMatGetData_Dense(A, f_array)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr)             :: A
    real(c_double), pointer :: f_array(:,:)

    ! internal variables
    type(c_ptr)     :: c_array
    integer(c_long) :: M, N

    !======= Internals ============

    ! get data pointer from N_Vector
    c_array = FSUNDenseMatrix_Data(A)

    ! get matrix size
    M = FSUNDenseMatrix_Rows(A)
    N = FSUNDenseMatrix_Columns(A)
    
    ! convert and reshape 1D data array
    call c_f_pointer(c_array, f_array, (/M,N/))

  end subroutine FSUNMatGetData_Dense

end module fsunmatrix_dense_mod
