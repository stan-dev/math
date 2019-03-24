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
! the SUNDIALS SPBCGS linear solver using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_spbcgs_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_SPBCGS(y, pretype, maxl) &
          bind(C,name='SUNLinSol_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FSUNLinSol_SPBCGS

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_SPBCGS(LS) &
          bind(C,name='SUNLinSolFree_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_SPBCGS

     ! =================================================================
     ! Setters
     ! =================================================================
     
     integer(c_int) function FSUNLinSol_SPBCGSSetPrecType(LS, pretype) &
          bind(C,name='SUNLinSol_SPBCGSSetPrecType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: pretype
     end function FSUNLinSol_SPBCGSSetPrecType
    
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

     integer(c_int) function FSUNLinSolGetType_SPBCGS(LS) &
          bind(C,name='SUNLinSolGetType_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_SPBCGS

     integer(c_int) function FSUNLinSolInitialize_SPBCGS(LS) &
          bind(C,name='SUNLinSolInitialize_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_SPBCGS
     
     integer(c_int) function FSUNLinSolSetATimes_SPBCGS(LS, A_data, ATimes) &
          bind(C,name='SUNLinSolSetATimes_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A_data
       type(c_funptr), value :: ATimes
     end function FSUNLinSolSetATimes_SPBCGS

     integer(c_int) function FSUNLinSolSetPreconditioner_SPBCGS(LS,     &
                                                             P_data, &
                                                             Pset,   &
                                                             Psol)   &
          bind(C,name='SUNLinSolSetPreconditioner_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: P_data
       type(c_funptr), value :: Pset
       type(c_funptr), value :: Psol
     end function FSUNLinSolSetPreconditioner_SPBCGS

     integer(c_int) function FSUNLinSolSetScalingVectors_SPBCGS(LS, s1, s2) &
          bind(C,name='SUNLinSolSetScalingVectors_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: s1
       type(c_ptr),    value :: s2
     end function FSUNLinSolSetScalingVectors_SPBCGS

     integer(c_int) function FSUNLinSolSetup_SPBCGS(LS, A) &
          bind(C,name='SUNLinSolSetup_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FSUNLinSolSetup_SPBCGS

     integer(c_int) function FSUNLinSolSolve_SPBCGS(LS, A, x, b, tol) &
          bind(C,name='SUNLinSolSolve_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_SPBCGS

     integer(c_int) function FSUNLinSolNumIters_SPBCGS(LS) &
          bind(C,name='SUNLinSolNumIters_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolNumIters_SPBCGS

     real(c_double) function FSUNLinSolResNorm_SPBCGS(LS) &
          bind(C,name='SUNLinSolResNorm_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResNorm_SPBCGS
     
     type(c_ptr) function FSUNLinSolResid_SPBCGS(LS) &
          bind(C,name='SUNLinSolResid_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResid_SPBCGS
     
     integer(c_long) function FSUNLinSolLastFlag_SPBCGS(LS) &
          bind(C,name='SUNLinSolLastFlag_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_SPBCGS

     integer(c_int) function FSUNLinSolSpace_SPBCGS(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_SPBCGS')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_SPBCGS

  end interface

end module fsunlinsol_spbcgs_mod
