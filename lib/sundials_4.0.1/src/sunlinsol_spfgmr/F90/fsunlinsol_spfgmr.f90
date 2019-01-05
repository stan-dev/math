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
! the SUNDIALS SPGMR linear solver using the ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunlinsol_spfgmr_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNLinSol_SPFGMR(y, pretype, maxl) &
          bind(C,name='SUNLinSol_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
       integer(c_int), value :: pretype
       integer(c_int), value :: maxl
     end function FSUNLinSol_SPFGMR

     ! =================================================================
     ! Destructors
     ! =================================================================

     subroutine FSUNLinSolFree_SPFGMR(LS) &
          bind(C,name='SUNLinSolFree_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end subroutine FSUNLinSolFree_SPFGMR

     ! =================================================================
     ! Setters
     ! =================================================================
     
     integer(c_int) function FSUNLinSol_SPFGMRSetPrecType(LS, pretype) &
          bind(C,name='SUNLinSol_SPFGMRSetPrecType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: pretype
     end function FSUNLinSol_SPFGMRSetPrecType
    
     integer(c_int) function FSUNLinSol_SPFGMRSetGSType(LS, gstype) &
          bind(C,name='SUNLinSol_SPFGMRSetGSType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: gstype
     end function FSUNLinSol_SPFGMRSetGSType

     integer(c_int) function FSUNLinSol_SPFGMRSetMaxRestarts(LS, maxrs) &
          bind(C,name='SUNLinSol_SPFGMRSetMaxRestarts')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       integer(c_int), value :: maxrs
     end function FSUNLinSol_SPFGMRSetMaxRestarts

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNLinSolGetType_SPFGMR(LS) &
          bind(C,name='SUNLinSolGetType_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolGetType_SPFGMR

     integer(c_int) function FSUNLinSolInitialize_SPFGMR(LS) &
          bind(C,name='SUNLinSolInitialize_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolInitialize_SPFGMR
     
     integer(c_int) function FSUNLinSolSetATimes_SPFGMR(LS, A_data, ATimes) &
          bind(C,name='SUNLinSolSetATimes_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A_data
       type(c_funptr), value :: ATimes
     end function FSUNLinSolSetATimes_SPFGMR

     integer(c_int) function FSUNLinSolSetPreconditioner_SPFGMR(LS,     &
                                                               P_data, &
                                                               Pset,   &
                                                               Psol)   &
          bind(C,name='SUNLinSolSetPreconditioner_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: P_data
       type(c_funptr), value :: Pset
       type(c_funptr), value :: Psol
     end function FSUNLinSolSetPreconditioner_SPFGMR

     integer(c_int) function FSUNLinSolSetScalingVectors_SPFGMR(LS, s1, s2) &
          bind(C,name='SUNLinSolSetScalingVectors_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: s1
       type(c_ptr),    value :: s2
     end function FSUNLinSolSetScalingVectors_SPFGMR

     integer(c_int) function FSUNLinSolSetup_SPFGMR(LS, A) &
          bind(C,name='SUNLinSolSetup_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FSUNLinSolSetup_SPFGMR

     integer(c_int) function FSUNLinSolSolve_SPFGMR(LS, A, x, b, tol) &
          bind(C,name='SUNLinSolSolve_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: LS
       type(c_ptr),    value :: A
       type(c_ptr),    value :: x
       type(c_ptr),    value :: b
       real(c_double), value :: tol
     end function FSUNLinSolSolve_SPFGMR

     integer(c_int) function FSUNLinSolNumIters_SPFGMR(LS) &
          bind(C,name='SUNLinSolNumIters_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolNumIters_SPFGMR

     real(c_double) function FSUNLinSolResNorm_SPFGMR(LS) &
          bind(C,name='SUNLinSolResNorm_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResNorm_SPFGMR
     
     type(c_ptr) function FSUNLinSolResid_SPFGMR(LS) &
          bind(C,name='SUNLinSolResid_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolResid_SPFGMR
     
     integer(c_long) function FSUNLinSolLastFlag_SPFGMR(LS) &
          bind(C,name='SUNLinSolLastFlag_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
     end function FSUNLinSolLastFlag_SPFGMR

     integer(c_int) function FSUNLinSolSpace_SPFGMR(LS, lenrwLS, leniwLS) &
          bind(C,name='SUNLinSolSpace_SPFGMR')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: LS
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FSUNLinSolSpace_SPFGMR

  end interface

end module fsunlinsol_spfgmr_mod
