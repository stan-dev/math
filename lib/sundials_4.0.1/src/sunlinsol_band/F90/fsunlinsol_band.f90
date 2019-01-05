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
! the SUNDIALS banded matrix using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_band_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_Band(y, A) &
          bind(C,name='SUNLinSol_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: y
       type(c_ptr), value :: A
     end function FSUNLinSol_Band

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_Band(LS) &
          bind(C,name='SUNLinSolFree_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_Band

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNLinSolGetType_Band(LS) &
          bind(C,name='SUNLinSolGetType_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_Band

     integer(c_int) function FSUNLinSolInitialize_Band(LS) &
          bind(C,name='SUNLinSolInitialize_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_Band

     integer(c_int) function FSUNLinSolSetup_Band(LS, A) &
          bind(C,name='SUNLinSolSetup_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FSUNLinSolSetup_Band

     integer(c_int) function FSUNLinSolSolve_Band(LS, A, x, b, tol) &
          bind(C,name='SUNLinSolSolve_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_Band

     integer(c_long) function FSUNLinSolLastFlag_Band(LS) &
          bind(C,name='SUNLinSolLastFlag_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_Band

     integer(c_int) function FSUNLinSolSpace_Band(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_Band')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_Band

  end interface

end module fsunlinsol_band_mod
