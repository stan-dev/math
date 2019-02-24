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
! the SUNDIALS sparse matrix using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunmatrix_sparse_mod

  use, intrinsic :: iso_c_binding, only : c_int

  integer(c_int), parameter :: CSC_MAT = 0
  integer(c_int), parameter :: CSR_MAT = 1

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNSparseMatrix(M, N, NNZ, sparsetype) &
          bind(C,name='SUNSparseMatrix')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_long), value :: M
       integer(c_long), value :: N
       integer(c_long), value :: NNZ
       integer(c_int),  value :: sparsetype
     end function FSUNSparseMatrix

     type(c_ptr) function FSUNSparseFromDenseMatrix(A, droptol, sparsetype) &
          bind(C,name='SUNSparseFromDenseMatrix')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value    :: A
       real(c_double), value :: droptol
       integer(c_int), value :: sparsetype
     end function FSUNSparseFromDenseMatrix

     type(c_ptr) function FSUNSparseFromBandMatrix(A, droptol, sparsetype) &
          bind(C,name='SUNSparseFromBandMatrix')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value    :: A
       real(c_double), value :: droptol
       integer(c_int), value :: sparsetype
     end function FSUNSparseFromBandMatrix

     integer(c_int) function FSUNSparseMatrix_Realloc(A) &
          bind(C,name='SUNSparseMatrix_Realloc')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_Realloc

     integer(c_int) function FSUNSparseMatrix_Reallocate(A, NNZ) &
          bind(C,name='SUNSparseMatrix_Reallocate')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value     :: A
       integer(c_long), value :: NNZ
     end function FSUNSparseMatrix_Reallocate

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNMatDestroy_Sparse(A) &
          bind(C,name='SUNMatDestroy_Sparse')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end subroutine FSUNMatDestroy_Sparse

     ! =================================================================
     ! Other routines
     ! =================================================================

     ! -----------------------------------------------------------------
     ! NOT INTERFACED SUNSparseMatrix_Print
     ! -----------------------------------------------------------------

     integer(c_long) function FSUNSparseMatrix_Rows(A) &
          bind(C,name='SUNSparseMatrix_Rows')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_Rows

     integer(c_long) function FSUNSparseMatrix_Columns(A) &
          bind(C,name='SUNSparseMatrix_Columns')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_Columns

     integer(c_long) function FSUNSparseMatrix_NNZ(A) &
          bind(C,name='SUNSparseMatrix_NNZ')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_NNZ

     integer(c_long) function FSUNSparseMatrix_NP(A) &
          bind(C,name='SUNSparseMatrix_NP')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_NP

     integer(c_int) function FSUNSparseMatrix_SparseType(A) &
          bind(C,name='SUNSparseMatrix_SparseType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_SparseType

     type(c_ptr) function FSUNSparseMatrix_Data(A) &
          bind(C,name='SUNSparseMatrix_Data')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_Data

     type(c_ptr) function FSUNSparseMatrix_IndexValues(A) &
          bind(C,name='SUNSparseMatrix_IndexValues')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_IndexValues

     type(c_ptr) function FSUNSparseMatrix_IndexPointers(A) &
          bind(C,name='SUNSparseMatrix_IndexPointers')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNSparseMatrix_IndexPointers

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNMatGetID_Sparse(A) &
          bind(C,name='SUNMatGetID_Sparse')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNMatGetID_Sparse

     type(c_ptr) function FSUNMatClone_Sparse(A) &
          bind(C,name='SUNMatClone_Sparse')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: A
     end function FSUNMatClone_Sparse

     integer(c_int) function FSUNMatZero_Sparse(A) &
          bind(C,name='SUNMatZero_Sparse')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
     end function FSUNMatZero_Sparse

     integer(c_int) function FSUNMatCopy_Sparse(A, B) &
          bind(C,name='SUNMatCopy_Sparse')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       type(c_ptr), value :: B
     end function FSUNMatCopy_Sparse

     integer(c_int) function FSUNMatScaleAdd_Sparse(c, A, B) &
          bind(C,name='SUNMatScaleAdd_Sparse')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: c
       type(c_ptr),    value :: A
       type(c_ptr),    value :: B
     end function FSUNMatScaleAdd_Sparse

     integer(c_int) function FSUNMatScaleAddI_Sparse(c, A) &
          bind(C,name='SUNMatScaleAddI_Sparse')
       use, intrinsic :: iso_c_binding
       real(c_double), value :: c
       type(c_ptr),    value :: A
     end function FSUNMatScaleAddI_Sparse

     integer(c_int) function FSUNMatMatvec_Sparse(A, x, y) &
          bind(C,name='SUNMatMatvec_Sparse')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       type(c_ptr), value :: x
       type(c_ptr), value :: y
     end function FSUNMatMatvec_Sparse

     integer(c_int) function FSUNMatSpace_Sparse(A, lenrw, leniw) &
          bind(C,name='SUNMatSpace_Sparse')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: A
       integer(c_long)    :: lenrw
       integer(c_long)    :: leniw
     end function FSUNMatSpace_Sparse

  end interface

contains

  ! ================================================================
  ! Helpful routines
  ! ================================================================

  ! ----------------------------------------------------------------
  ! FSUNMatGetData_Sparse
  !
  ! Extracts data array from a SUNDIALS Sparse Matrix
  ! ----------------------------------------------------------------

  subroutine FSUNMatGetData_Sparse(A, f_data, f_idxval, f_idxptr)

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none

    ! calling variables
    type(c_ptr)              :: A
    real(c_double),  pointer :: f_data(:)
    integer(c_long), pointer :: f_idxval(:)
    integer(c_long), pointer :: f_idxptr(:)

    ! internal variables
    type(c_ptr) :: c_data
    type(c_ptr) :: c_idxval
    type(c_ptr) :: c_idxptr

    integer(c_long) :: NNZ
    integer(c_long) :: NP

    !======= Internals ============

    ! get data pointer from N_Vector
    c_data   = FSUNSparseMatrix_Data(A)
    c_idxval = FSUNSparseMatrix_IndexValues(A)
    c_idxptr = FSUNSparseMatrix_IndexPointers(A)

    ! get matrix size
    NNZ = FSUNSparseMatrix_NNZ(A)
    NP  = FSUNSparseMatrix_NP(A)

    ! convert and reshape 1D data array
    call c_f_pointer(c_data,   f_data,   (/NNZ/))
    call c_f_pointer(c_idxval, f_idxval, (/NNZ/))
    call c_f_pointer(c_idxptr, f_idxptr, (/NP+1/))

  end subroutine FSUNMatGetData_Sparse

end module fsunmatrix_sparse_mod
