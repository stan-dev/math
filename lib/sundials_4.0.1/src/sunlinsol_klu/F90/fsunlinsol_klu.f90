! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! LLNS Copyright Start
! Copyright (c) 2017, Lawrence Livermore National Security
! This work was performed under the auspices of the U.S. Department
! of Energy by Lawrence Livermore National Laboratory in part under
! Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
! Produced at the Lawrence Livermore National Laboratory.
! All rights reserved.
! For details, see the LICENSE file.
! LLNS Copyright End
! -----------------------------------------------------------------
! This file contains a Fortran module for interfacing directly with
! the SUNDIALS sparse matrix using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_klu_mod

  use, intrinsic :: iso_c_binding, only : c_int
  
  integer(c_int), parameter :: SUNKLU_ORDERING_DEFAULT = 1 ! COLAMD
  integer(c_int), parameter :: SUNKLU_REINIT_FULL      = 1
  integer(c_int), parameter :: SUNKLU_REINIT_PARTIAL   = 2

  !======= Interfaces ========
  interface

    ! =================================================================
    ! Constructors
    ! =================================================================

    type(c_ptr) function FSUNKLU(y, A) &
        bind(C,name='SUNKLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: y
      type(c_ptr), value :: A
    end function FSUNKLU

    ! =================================================================
    ! Destructors
    ! =================================================================

    subroutine FSUNLinSolFree_KLU(LS) &
        bind(C,name='SUNLinSolFree_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end subroutine FSUNLinSolFree_KLU

    ! =================================================================
    ! Setter/init routines
    ! =================================================================

    integer(c_int) function FSUNKLUReInit(LS, A, nnz, reinit_type) &
      bind(C,name='SUNKLUReInit')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value     :: LS
      type(c_ptr), value     :: A
      integer(c_long), value :: nnz
      integer(c_int), value  :: reinit_type
    end function FSUNKLUReInit

    integer(c_int) function FSUNKLUSetOrdering(LS, ordering_choice) &
        bind(C,name='SUNKLUSetOrdering')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: LS
      integer(c_int), value :: ordering_choice
    end function FSUNKLUSetOrdering

    ! =================================================================
    ! Operations
    ! =================================================================

    integer(c_int) function FSUNLinSolGetType_KLU(LS) &
        bind(C,name='SUNLinSolGetType_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end function FSUNLinSolGetType_KLU

    integer(c_int) function FSUNLinSolInitialize_KLU(LS) &
        bind(C,name='SUNLinSolInitialize_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end function FSUNLinSolInitialize_KLU

    integer(c_int) function FSUNLinSolSetup_KLU(LS, A) &
        bind(C,name='SUNLinSolSetup_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
      type(c_ptr), value :: A
    end function FSUNLinSolSetup_KLU

    integer(c_int) function FSUNLinSolSolve_KLU(LS, A, x, b, tol) &
        bind(C,name='SUNLinSolSolve_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: LS
      type(c_ptr), value    :: A
      type(c_ptr), value    :: x
      type(c_ptr), value    :: b
      real(c_double), value :: tol
    end function FSUNLinSolSolve_KLU

    integer(c_long) function FSUNLinSolLastFlag_KLU(LS) &
        bind(C,name='SUNLinSolLastFlag_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
    end function FSUNLinSolLastFlag_KLU

    integer(c_int) function FSUNLinSolSpace_KLU(LS, lenrwLS, leniwLS) &
        bind(C,name='SUNLinSolSpace_KLU')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: LS
      integer(c_long), value :: lenrwLS
      integer(c_long), value :: leniwLS
    end function FSUNLinSolSpace_KLU

  end interface

end module fsunlinsol_klu_mod
