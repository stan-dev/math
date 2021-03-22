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
! This is an example custom SUNLINEARSOLVER module for reaction
! problems.  Since these include no inter-cell coupling, the
! Jacobian is block-diagonal, with N blocks of size (Nvar x Nvar)
! each.  The fsunmatrix_fortran_mod module stores these blocks in a
! 3-dimensional Fortran array.  This structure is leveraged here for
! an efficient linear solve.
! ------------------------------------------------------------------

module fsunlinsol_fortran_mod

  use, intrinsic :: iso_c_binding
  use fsundials_linearsolver_mod
  use fsunmatrix_fortran_mod
  use fnvector_fortran_mod

  implicit none

  ! ----------------------------------------------------------------
  type, public :: FLinSol
     integer(c_long)              :: Nvar
     integer(c_long)              :: N
     integer(c_long), allocatable :: pivots(:,:)
  end type FLinSol
  ! ----------------------------------------------------------------

contains

  ! ----------------------------------------------------------------
  function FSUNLinSolNew_Fortran(Nvar, N) result(sunls_S)

    implicit none
    integer(c_long),           value   :: Nvar
    integer(c_long),           value   :: N
    type(SUNLinearSolver),     pointer :: sunls_S
    type(SUNLinearSolver_Ops), pointer :: ops
    type(FLinSol),             pointer :: content

    ! allocate output SUNLinearSolver structure
    sunls_S => FSUNLinSolNewEmpty()

    ! allocate and fill content structure
    allocate(content)
    allocate(content%pivots(Nvar,N))
    content%Nvar = NVar
    content%N    = N

    ! attach the content structure to the output SUNMatrix
    sunls_S%content = c_loc(content)

    ! access the ops structure, and set internal function pointers
    call c_f_pointer(sunls_S%ops, ops)
    ops%gettype = c_funloc(FSUNLinSolGetType_Fortran)
    ops%setup   = c_funloc(FSUNLinSolSetup_Fortran)
    ops%solve   = c_funloc(FSUNLinSolSolve_Fortran)
    ops%space   = c_funloc(FSUNLinSolSpace_Fortran)
    ops%free    = c_funloc(FSUNLinSolFree_Fortran)

  end function FSUNLinSolNew_Fortran

  ! ----------------------------------------------------------------
  function FSUNLinSolGetFLinSol(sunls_S) result(S)

    implicit none
    type(SUNLinearSolver)  :: sunls_S
    type(FLinSol), pointer :: S

    ! extract Fortran matrix structure to output
    call c_f_pointer(sunls_S%content, S)

    return

  end function FSUNLinSolGetFLinSol

  ! ----------------------------------------------------------------
  integer(SUNLinearSolver_Type) function FSUNLinSolGetType_Fortran(sunls_S) &
    result(id) bind(C)

    implicit none
    type(SUNLinearSolver) :: sunls_S

    id = SUNLINEARSOLVER_DIRECT
    return

  end function FSUNLinSolGetType_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNLinSolFree_Fortran(sunls_S) &
       result(ierr) bind(C)

    implicit none
    type(SUNLinearSolver), target  :: sunls_S
    type(FLinSol),         pointer :: S

    ! access FLinSol structure
    S => FSUNLinSolGetFLinSol(sunls_S)

    ! deallocate pivots
    deallocate(S%pivots)

    ! deallocate the underlying Fortran object (the content)
    deallocate(S)

    ! set SUNLinearSolver structure members to NULL and return
    sunls_S%content = C_NULL_PTR

    ! deallocate overall SUNLinearSolver structure
    call FSUNLinSolFreeEmpty(sunls_S)

    ! return with success
    ierr = 0
    return

  end function FSUNLinSolFree_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNLinSolSetup_Fortran(sunls_S, sunmat_A) &
       result(ierr) bind(C)

    implicit none
    type(SUNLinearSolver)   :: sunls_S
    type(SUNMatrix)         :: sunmat_A
    type(FLinSol),  pointer :: S
    type(FMat),     pointer :: AMat
    integer(c_long)         :: i, j, k, l
    real(c_double)          :: temp
    real(c_double), pointer :: A(:,:)

    ! extract Fortran structures to work with
    S => FSUNLinSolGetFLinSol(sunls_S)
    AMat => FSUNMatGetFMat(sunmat_A)

    ! perform LU factorization of each block on diagonal
    do i = 1,S%N

       ! set 2D pointer to this diagonal block
       A => AMat%data(:,:,i)

       ! k-th elimination step number
       do k = 1,S%Nvar

          ! find l = pivot row number
          l = k
          do j = k+1,S%Nvar
             if (dabs(A(j,k)) > dabs(A(l,k))) then
                l = j
             end if
          end do
          S%pivots(k,i) = l

          ! check for zero pivot element
          if (A(l,k) == 0.d0) then
             ierr = int(k, c_int)
             return
          end if

          ! swap a(k,1:n) and a(l,1:n) if necessary
          if ( l /= k ) then
             do j = 1,S%Nvar
                temp = A(l,j)
                A(l,j) = A(k,j)
                A(k,j) = temp
             end do
          end if

          ! Scale the elements below the diagonal in
          ! column k by 1.0/a(k,k). After the above swap
          ! a(k,k) holds the pivot element. This scaling
          ! stores the pivot row multipliers a(i,k)/a(k,k)
          ! in a(i,k), i=k+1, ..., m-1.
          A(k+1:S%Nvar,k) = A(k+1:S%Nvar,k) / A(k,k)

          ! row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., Nvar
          ! row k is the pivot row after swapping with row l.
          ! The computation is done one column at a time
          do j = k+1,S%Nvar
             if (A(k,j) /= 0.d0) then
                A(k+1:S%Nvar,j) = A(k+1:S%Nvar,j) - A(k,j) * A(k+1:S%Nvar,k)
             end if
          end do

       end do
    end do

    ! return with success
    ierr = 0
    return

  end function FSUNLinSolSetup_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNLinSolSolve_Fortran(sunls_S, sunmat_A, &
       sunvec_x, sunvec_b, tol) result(ierr) bind(C)

    implicit none
    type(SUNLinearSolver)   :: sunls_S
    type(SUNMatrix)         :: sunmat_A
    type(N_Vector)          :: sunvec_x
    type(N_Vector)          :: sunvec_b
    real(c_double), value   :: tol
    type(FLinSol),  pointer :: S
    type(FMat),     pointer :: AMat
    type(FVec),     pointer :: xvec, bvec
    integer(c_long)         :: i, k, pk
    real(c_double)          :: temp
    real(c_double), pointer :: A(:,:), x(:)

    ! extract Fortran structures to work with
    S => FSUNLinSolGetFLinSol(sunls_S)
    AMat => FSUNMatGetFMat(sunmat_A)
    xvec => FN_VGetFVec(sunvec_x)
    bvec => FN_VGetFVec(sunvec_b)

    ! copy b into x
    xvec%data(:,:) = bvec%data(:,:)

    ! perform solve using LU-factored blocks on matrix diagonal
    do i = 1,S%N

       ! set pointer to this block of overall linear system
       A => AMat%data(:,:,i)
       x => xvec%data(:,i)

       ! Permute x, based on pivot information in p
       do k = 1,S%Nvar
          pk = S%pivots(k,i)
          if (pk /= k) then
             temp = x(k)
             x(k) = x(pk)
             x(pk) = temp
          end if
       end do

       ! Solve Ly = x, store solution y in x
       do k = 1,S%Nvar-1
          x(k+1:S%Nvar) = x(k+1:S%Nvar) - x(k)*A(k+1:S%Nvar,k)
       end do

       ! Solve Ux = y (y is initially stored in x)
       do k = S%Nvar,2,-1
          x(k) = x(k)/A(k,k)
          x(1:k-1) = x(1:k-1) - A(1:k-1,k)*x(k)
       end do
       x(1) = x(1)/A(1,1)

    end do

    ! return with success
    ierr = 0
    return

  end function FSUNLinSolSolve_Fortran

  ! ----------------------------------------------------------------
  integer(c_int) function FSUNLinSolSpace_Fortran(sunls_S, lrw, liw) &
       result(ierr) bind(C)

    implicit none
    type(SUNLinearSolver)  :: sunls_S
    integer(c_int64_t)     :: lrw(1)
    integer(c_int64_t)     :: liw(1)
    type(FLinSol), pointer :: S

    ! extract Fortran structure to work with
    S => FSUNLinSolGetFLinSol(sunls_S)

    ! set output arguments and return
    lrw(1) = (S%Nvar)*(S%N)
    liw(1) = 2

    ! return with success
    ierr = 0
    return

  end function FSUNLinSolSpace_Fortran

end module fsunlinsol_fortran_mod
! ------------------------------------------------------------------
