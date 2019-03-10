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
! the SUNDIALS SPTFQMR linear solver using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_sptfqmr_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_SPTFQMR(y, pretype, maxl) &
          bind(C,name='SUNLinSol_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FSUNLinSol_SPTFQMR

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_SPTFQMR(LS) &
          bind(C,name='SUNLinSolFree_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_SPTFQMR

     ! =================================================================
     ! Setters
     ! =================================================================
     
     integer(c_int) function FSUNLinSol_SPTFQMRSetPrecType(LS, pretype) &
          bind(C,name='SUNLinSol_SPTFQMRSetPrecType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: pretype
     end function FSUNLinSol_SPTFQMRSetPrecType

     integer(c_int) function FSUNLinSol_SPTFQMRSetMaxl(LS, maxl) &
          bind(C,name='SUNLinSol_SPTFQMRSetMaxl')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: maxl
     end function FSUNLinSol_SPTFQMRSetMaxl

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNLinSolGetType_SPTFQMR(LS) &
          bind(C,name='SUNLinSolGetType_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_SPTFQMR

     integer(c_int) function FSUNLinSolInitialize_SPTFQMR(LS) &
          bind(C,name='SUNLinSolInitialize_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_SPTFQMR
     
     integer(c_int) function FSUNLinSolSetATimes_SPTFQMR(LS, A_data, ATimes) &
          bind(C,name='SUNLinSolSetATimes_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A_data
       type(c_funptr), value :: ATimes
     end function FSUNLinSolSetATimes_SPTFQMR

     integer(c_int) function FSUNLinSolSetPreconditioner_SPTFQMR(LS,     &
                                                               P_data, &
                                                               Pset,   &
                                                               Psol)   &
          bind(C,name='SUNLinSolSetPreconditioner_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: P_data
       type(c_funptr), value :: Pset
       type(c_funptr), value :: Psol
     end function FSUNLinSolSetPreconditioner_SPTFQMR

     integer(c_int) function FSUNLinSolSetScalingVectors_SPTFQMR(LS, s1, s2) &
          bind(C,name='SUNLinSolSetScalingVectors_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: s1
       type(c_ptr),    value :: s2
     end function FSUNLinSolSetScalingVectors_SPTFQMR

     integer(c_int) function FSUNLinSolSetup_SPTFQMR(LS, A) &
          bind(C,name='SUNLinSolSetup_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FSUNLinSolSetup_SPTFQMR

     integer(c_int) function FSUNLinSolSolve_SPTFQMR(LS, A, x, b, tol) &
          bind(C,name='SUNLinSolSolve_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_SPTFQMR

     integer(c_int) function FSUNLinSolNumIters_SPTFQMR(LS) &
          bind(C,name='SUNLinSolNumIters_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolNumIters_SPTFQMR

     real(c_double) function FSUNLinSolResNorm_SPTFQMR(LS) &
          bind(C,name='SUNLinSolResNorm_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResNorm_SPTFQMR
     
     type(c_ptr) function FSUNLinSolResid_SPTFQMR(LS) &
          bind(C,name='SUNLinSolResid_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResid_SPTFQMR
     
     integer(c_long) function FSUNLinSolLastFlag_SPTFQMR(LS) &
          bind(C,name='SUNLinSolLastFlag_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_SPTFQMR

     integer(c_int) function FSUNLinSolSpace_SPTFQMR(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_SPTFQMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_SPTFQMR

  end interface

end module fsunlinsol_sptfqmr_mod
