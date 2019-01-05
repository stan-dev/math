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
! the SUNDIALS dense linear solver using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_dense_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_Dense(y, A) &
          bind(C,name='SUNLinSol_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: y
       type(c_ptr), value :: A
     end function FSUNLinSol_Dense
     
     ! Deprecated
     type(c_ptr) function FSUNDenseLinearSolver(y, A) &
          bind(C,name='SUNDenseLinearSolver')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: y
       type(c_ptr), value :: A
     end function FSUNDenseLinearSolver

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_Dense(LS) &
          bind(C,name='SUNLinSolFree_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_Dense

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNLinSolGetType_Dense(LS) &
          bind(C,name='SUNLinSolGetType_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_Dense

     integer(c_int) function FSUNLinSolInitialize_Dense(LS) &
          bind(C,name='SUNLinSolInitialize_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_Dense

     integer(c_int) function FSUNLinSolSetup_Dense(LS, A) &
          bind(C,name='SUNLinSolSetup_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FSUNLinSolSetup_Dense

     integer(c_int) function FSUNLinSolSolve_Dense(LS, A, x, b, tol) &
          bind(C,name='SUNLinSolSolve_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_Dense

     integer(c_long) function FSUNLinSolLastFlag_Dense(LS) &
          bind(C,name='SUNLinSolLastFlag_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_Dense

     integer(c_int) function FSUNLinSolSpace_Dense(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_Dense')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_Dense

  end interface

end module fsunlinsol_dense_mod
