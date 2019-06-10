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
! the SUNDIALS SPGMR linear solver using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_spgmr_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_SPGMR(y, pretype, maxl) &
          bind(C,name='SUNLinSol_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FSUNLinSol_SPGMR

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_SPGMR(LS) &
          bind(C,name='SUNLinSolFree_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_SPGMR

     ! =================================================================
     ! Setters
     ! =================================================================
     
     integer(c_int) function FSUNLinSol_SPGMRSetPrecType(LS, pretype) &
          bind(C,name='SUNLinSol_SPGMRSetPrecType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: pretype
     end function FSUNLinSol_SPGMRSetPrecType
    
     integer(c_int) function FSUNLinSol_SPGMRSetGSType(LS, gstype) &
          bind(C,name='SUNLinSol_SPGMRSetGSType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: gstype
     end function FSUNLinSol_SPGMRSetGSType

     integer(c_int) function FSUNLinSol_SPGMRSetMaxRestarts(LS, maxrs) &
          bind(C,name='SUNLinSol_SPGMRSetMaxRestarts')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: maxrs
     end function FSUNLinSol_SPGMRSetMaxRestarts

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNLinSolGetType_SPGMR(LS) &
          bind(C,name='SUNLinSolGetType_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_SPGMR

     integer(c_int) function FSUNLinSolInitialize_SPGMR(LS) &
          bind(C,name='SUNLinSolInitialize_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_SPGMR
     
     integer(c_int) function FSUNLinSolSetATimes_SPGMR(LS, A_data, ATimes) &
          bind(C,name='SUNLinSolSetATimes_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A_data
       type(c_funptr), value :: ATimes
     end function FSUNLinSolSetATimes_SPGMR

     integer(c_int) function FSUNLinSolSetPreconditioner_SPGMR(LS,     &
                                                               P_data, &
                                                               Pset,   &
                                                               Psol)   &
          bind(C,name='SUNLinSolSetPreconditioner_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: P_data
       type(c_funptr), value :: Pset
       type(c_funptr), value :: Psol
     end function FSUNLinSolSetPreconditioner_SPGMR

     integer(c_int) function FSUNLinSolSetScalingVectors_SPGMR(LS, s1, s2) &
          bind(C,name='SUNLinSolSetScalingVectors_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: s1
       type(c_ptr),    value :: s2
     end function FSUNLinSolSetScalingVectors_SPGMR

     integer(c_int) function FSUNLinSolSetup_SPGMR(LS, A) &
          bind(C,name='SUNLinSolSetup_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FSUNLinSolSetup_SPGMR

     integer(c_int) function FSUNLinSolSolve_SPGMR(LS, A, x, b, tol) &
          bind(C,name='SUNLinSolSolve_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_SPGMR

     integer(c_int) function FSUNLinSolNumIters_SPGMR(LS) &
          bind(C,name='SUNLinSolNumIters_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolNumIters_SPGMR

     real(c_double) function FSUNLinSolResNorm_SPGMR(LS) &
          bind(C,name='SUNLinSolResNorm_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResNorm_SPGMR
     
     type(c_ptr) function FSUNLinSolResid_SPGMR(LS) &
          bind(C,name='SUNLinSolResid_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResid_SPGMR
     
     integer(c_long) function FSUNLinSolLastFlag_SPGMR(LS) &
          bind(C,name='SUNLinSolLastFlag_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_SPGMR

     integer(c_int) function FSUNLinSolSpace_SPGMR(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_SPGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_SPGMR

  end interface

end module fsunlinsol_spgmr_mod
