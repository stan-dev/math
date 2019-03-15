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
! the SUNDIALS PCG linear solver using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_pcg_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_PCG(y, pretype, maxl) &
          bind(C,name='SUNLinSol_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FSUNLinSol_PCG

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_PCG(LS) &
          bind(C,name='SUNLinSolFree_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_PCG

     ! =================================================================
     ! Setters
     ! =================================================================
     
     integer(c_int) function FSUNLinSol_PCGSetPrecType(LS, pretype) &
          bind(C,name='SUNLinSol_PCGSetPrecType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: pretype
     end function FSUNLinSol_PCGSetPrecType
    
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

     integer(c_int) function FSUNLinSolGetType_PCG(LS) &
          bind(C,name='SUNLinSolGetType_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_PCG

     integer(c_int) function FSUNLinSolInitialize_PCG(LS) &
          bind(C,name='SUNLinSolInitialize_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_PCG
     
     integer(c_int) function FSUNLinSolSetATimes_PCG(LS, A_data, ATimes) &
          bind(C,name='SUNLinSolSetATimes_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A_data
       type(c_funptr), value :: ATimes
     end function FSUNLinSolSetATimes_PCG

     integer(c_int) function FSUNLinSolSetPreconditioner_PCG(LS,     &
                                                             P_data, &
                                                             Pset,   &
                                                             Psol)   &
          bind(C,name='SUNLinSolSetPreconditioner_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: P_data
       type(c_funptr), value :: Pset
       type(c_funptr), value :: Psol
     end function FSUNLinSolSetPreconditioner_PCG

     integer(c_int) function FSUNLinSolSetScalingVectors_PCG(LS, s1, nul) &
          bind(C,name='SUNLinSolSetScalingVectors_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: s1
       type(c_ptr),    value :: nul
     end function FSUNLinSolSetScalingVectors_PCG

     integer(c_int) function FSUNLinSolSetup_PCG(LS, nul) &
          bind(C,name='SUNLinSolSetup_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: nul
     end function FSUNLinSolSetup_PCG

     integer(c_int) function FSUNLinSolSolve_PCG(LS, nul, x, b, tol) &
          bind(C,name='SUNLinSolSolve_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: nul
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_PCG

     integer(c_int) function FSUNLinSolNumIters_PCG(LS) &
          bind(C,name='SUNLinSolNumIters_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolNumIters_PCG

     real(c_double) function FSUNLinSolResNorm_PCG(LS) &
          bind(C,name='SUNLinSolResNorm_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResNorm_PCG
     
     type(c_ptr) function FSUNLinSolResid_PCG(LS) &
          bind(C,name='SUNLinSolResid_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResid_PCG
     
     integer(c_long) function FSUNLinSolLastFlag_PCG(LS) &
          bind(C,name='SUNLinSolLastFlag_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_PCG

     integer(c_int) function FSUNLinSolSpace_PCG(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_PCG')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_PCG

  end interface

end module fsunlinsol_pcg_mod
