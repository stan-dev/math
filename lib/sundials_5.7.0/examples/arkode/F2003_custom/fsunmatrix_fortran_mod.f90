! ------------------------------------------------------------------
! Programmer(s): Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! This is an example custom SUNMATRIX module for reaction problems.
! Since these include no inter-cell coupling, the Jacobian is
! block-diagonal, with N blocks of size (Nvar x Nvar) each.  We
! store these blocks in a 3-dimensional Fortran array.  This
! structure is subsequently leveraged in the custom
! fsunlinsol_fortran_mod module.
! ------------------------------------------------------------------

module fsunmatrix_fortran_mod

  use, intrinsic :: iso_c_binding
  use fsundials_matrix_mod
  use fnvector_fortran_mod

  implicit none

  ! ----------------------------------------------------------------
  type, public :: FMat
    logical                 :: own_data
    integer(c_long)         :: Nvar
    integer(c_long)         :: N
    real(c_double), pointer :: data(:,:,:)
  end type FMat
  ! ----------------------------------------------------------------

contains

  ! ----------------------------------------------------------------
  function FSUNMatNew_Fortran(Nvar, N) result(sunmat_A)

    implicit none
    integer(c_long),     value   :: Nvar
    integer(c_long),     value   :: N
    type(SUNMatrix),     pointer :: sunmat_A
    type(SUNMatrix_Ops), pointer :: ops
    type(FMat),          pointer :: content

    ! allocate output SUNMatrix structure
    sunmat_A => FSUNMatNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    allocate(content%data(Nvar,Nvar,N))
    content%own_data = .true.
    content%Nvar     = NVar
    content%N        = N

    ! attach the content structure to the output SUNMatrix
    sunmat_A%content = c_loc(content)

    ! access the SUNMatrix ops structure, and set internal function pointers
    call c_f_pointer(sunmat_A%ops, ops)
    ops%getid     = c_funloc(FSUNMatGetID_Fortran)
    ops%clone     = c_funloc(FSUNMatClone_Fortran)
    ops%destroy   = c_funloc(FSUNMatDestroy_Fortran)
    ops%zero      = c_funloc(FSUNMatZero_Fortran)
    ops%copy      = c_funloc(FSUNMatCopy_Fortran)
    ops%scaleadd  = c_funloc(FSUNMatScaleAdd_Fortran)
    ops%scaleaddi = c_funloc(FSUNMatScaleAddI_Fortran)
    ops%matvec    = c_funloc(FSUNMatMatvec_Fortran)
    ops%space     = c_funloc(FSUNMatSpace_Fortran)

  end function FSUNMatNew_Fortran

  ! ----------------------------------------------------------------
  function FSUNMatGetFMat(sunmat_A) result(A)

    implicit none
    type(SUNMatrix)     :: sunmat_A
    type(FMat), pointer :: A

    ! extract Fortran matrix structure to output
    call c_f_pointer(sunmat_A%content, A)

    return

  end function FSUNMatGetFMat

  ! ----------------------------------------------------------------
  integer(SUNMatrix_ID) function FSUNMatGetID_Fortran(sunmat_A) &
    result(id) bind(C)

    implicit none
    type(SUNMatrix) :: sunmat_A

    id = SUNMATRIX_CUSTOM
    return

  end function FSUNMatGetID_Fortran

  ! ----------------------------------------------------------------
  function FSUNMatClone_Fortran(sunmat_A) result(B_ptr) bind(C)

    implicit none
    type(SUNMatrix)          :: sunmat_A
    type(SUNMatrix), pointer :: sunmat_B
    type(c_ptr)              :: B_ptr
    integer(c_int)           :: retval
    type(FMat),      pointer :: A, B

    ! extract Fortran matrix structure to work with
    A => FSUNMatGetFMat(sunmat_A)

    ! allocate output N_Vector structure
    sunmat_B => FSUNMatNewEmpty()

    ! copy operations from x into y
    retval = FSUNMatCopyOps(sunmat_A, sunmat_B)

    ! allocate and clone content structure
    allocate(B)
    allocate(B%data(A%Nvar,A%Nvar,A%N))
    B%own_data = .true.
    B%Nvar = A%Nvar
    B%N = A%N

    ! attach the content structure to the output N_Vector
    sunmat_B%content = c_loc(B)

    ! set the c_ptr output
    B_ptr = c_loc(sunmat_B)
    return

  end function FSUNMatClone_Fortran

  ! ----------------------------------------------------------------
  subroutine FSUNMatDestroy_Fortran(sunmat_A) bind(C)

    implicit none
    type(SUNMatrix), target  :: sunmat_A
    type(FMat),      pointer :: A

    ! access FMat structure
    A => FSUNMatGetFMat(sunmat_A)

    ! if matrix owns the data, then deallocate
    if (A%own_data)  deallocate(A%data)

    ! deallocate the underlying Fortran object (the content)
    deallocate(A)

    ! set SUNMatrix structure members to NULL and return
    sunmat_A%content = C_NULL_PTR

    ! deallocate overall SUNMatrix structure
    call FSUNMatFreeEmpty(sunmat_A)

    return

  end subroutine FSUNMatDestroy_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNMatZero_Fortran(sunmat_A) &
       result(ierr) bind(C)

    implicit none
    type(SUNMatrix)     :: sunmat_A
    type(FMat), pointer :: A

    ! extract Fortran matrix structure to work with
    A => FSUNMatGetFMat(sunmat_A)

    ! set all entries to zero (whole array operation)
    A%data(:,:,:) = 0.d0

    ! return with success
    ierr = 0
    return

  end function FSUNMatZero_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNMatCopy_Fortran(sunmat_A, sunmat_B) &
       result(ierr) bind(C)

    implicit none
    type(SUNMatrix)     :: sunmat_A
    type(SUNMatrix)     :: sunmat_B
    type(FMat), pointer :: A, B

    ! extract Fortran matrix structures to work with
    A => FSUNMatGetFMat(sunmat_A)
    B => FSUNMatGetFMat(sunmat_B)

    ! copy all entries from A into B (whole array operation)
    B%data(:,:,:) = A%data(:,:,:)

    ! return with success
    ierr = 0
    return

  end function FSUNMatCopy_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNMatScaleAdd_Fortran(c, sunmat_A, sunmat_B) &
       result(ierr) bind(C)

    implicit none
    real(c_double), value   :: c
    type(SUNMatrix)         :: sunmat_A
    type(SUNMatrix)         :: sunmat_B
    type(FMat),     pointer :: A, B

    ! extract Fortran matrix structures to work with
    A => FSUNMatGetFMat(sunmat_A)
    B => FSUNMatGetFMat(sunmat_B)

    ! A = c*A + B (whole array operation)
    A%data(:,:,:) = c * A%data(:,:,:) + B%data(:,:,:)

    ! return with success
    ierr = 0
    return

  end function FSUNMatScaleAdd_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNMatScaleAddI_Fortran(c, sunmat_A) &
       result(ierr) bind(C)

    implicit none
    real(c_double), value :: c
    type(SUNMatrix)       :: sunmat_A
    type(FMat), pointer   :: A
    integer(c_long)       :: i, j, k

    ! extract Fortran matrix structure to work with
    A => FSUNMatGetFMat(sunmat_A)

    ! A = c*A + I
    do k = 1,A%N
       do j = 1,A%Nvar
          do i = 1,A%Nvar
             A%data(i,j,k) = c * A%data(i,j,k)
          end do
          A%data(j,j,k) = A%data(j,j,k) + 1.d0
       end do
    end do

    ! return with success
    ierr = 0
    return

  end function FSUNMatScaleAddI_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNMatMatvec_Fortran(sunmat_A, sunvec_x, sunvec_y) &
       result(ierr) bind(C)

    implicit none
    type(SUNMatrix)     :: sunmat_A
    type(N_Vector)      :: sunvec_x
    type(N_Vector)      :: sunvec_y
    type(FMat), pointer :: A
    type(FVec), pointer :: x, y
    integer(c_long)     :: i

    ! extract Fortran matrix and vector structures to work with
    A => FSUNMatGetFMat(sunmat_A)
    x => FN_VGetFVec(sunvec_x)
    y => FN_VGetFVec(sunvec_y)

    ! y = A*x
    do i = 1,A%N
       y%data(:,i) = matmul(A%data(:,:,i), x%data(:,i))
    end do

    ! return with success
    ierr = 0
    return

  end function FSUNMatMatvec_Fortran

  ! ----------------------------------------------------------------
  subroutine FSUNMatSpace_Fortran(sunmat_A, lrw, liw) bind(C)

    implicit none
    type(SUNMatrix)     :: sunmat_A
    integer(c_int64_t)  :: lrw(1)
    integer(c_int64_t)  :: liw(1)
    type(FMat), pointer :: A

    ! extract Fortran matrix structure to work with
    A => FSUNMatGetFMat(sunmat_A)

    ! set output arguments and return
    lrw(1) = (A%Nvar)*(A%Nvar)*(A%N)
    liw(1) = 3
    return

  end subroutine FSUNMatSpace_Fortran

end module fsunmatrix_fortran_mod
! ------------------------------------------------------------------
