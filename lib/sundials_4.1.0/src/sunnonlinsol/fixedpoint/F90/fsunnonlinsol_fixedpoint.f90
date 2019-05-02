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
! the SUNDIALS fixed-point nonlinear solver using the ISO_C_BINDING
! module.
! -----------------------------------------------------------------

module fsunnonlinsol_fixedpoint_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNNonlinSol_FixedPoint(y, m) &
         bind(C,name='SUNNonlinSol_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
       integer(c_int), value :: m
     end function FSUNNonlinSol_FixedPoint

     type(c_ptr) function FSUNNonlinSol_FixedPointSens(cnt, y, m) &
         bind(C,name='SUNNonlinSol_FixedPointSens')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: cnt
       type(c_ptr),    value :: y
       integer(c_int), value :: m
     end function FSUNNonlinSol_FixedPointSens

     ! =================================================================
     ! Destructors
     ! =================================================================

     integer(c_int) function FSUNNonlinSolFree_FixedPoint(NLS) &
         bind(C,name='SUNNonlinSolFree_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolFree_FixedPoint

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNNonlinSolGetType_FixedPoint(NLS) &
         bind(C,name='SUNNonlinSolGetType_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolGetType_FixedPoint

     integer(c_int) function FSUNNonlinSolInitialize_FixedPoint(NLS) &
         bind(C,name='SUNNonlinSolInitialize_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolInitialize_FixedPoint

     integer(c_int) function FSUNNonlinSolSolve_FixedPoint(NLS, y0, y, w, tol, &
                                                           callSetup, mem) &
         bind(C,name='SUNNonlinSolSolve_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_ptr),    value :: y0
       type(c_ptr),    value :: y
       type(c_ptr),    value :: w
       real(c_double), value :: tol
       integer(c_int), value :: callSetup
       type(c_ptr),    value :: mem
     end function FSUNNonlinSolSolve_FixedPoint
    
    ! =================================================================
    ! Set functions
    ! =================================================================

     integer(c_int) function FSUNNonlinSolSetSysFn_FixedPoint(NLS, SysFn) &
         bind(C,name='SUNNonlinSolSetSysFn_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: SysFn
     end function FSUNNonlinSolSetSysFn_FixedPoint

     integer(c_int) function FSUNNonlinSolSetConvTestFn_FixedPoint(NLS, CTestFN) &
         bind(C,name='SUNNonlinSolSetConvTestFn_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: CTestFN
     end function FSUNNonlinSolSetConvTestFn_FixedPoint

     integer(c_int) function FSUNNonlinSolSetMaxIters_FixedPoint(NLS, maxiters) &
         bind(C,name='SUNNonlinSolSetMaxIters_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       integer(c_int), value :: maxiters
     end function FSUNNonlinSolSetMaxIters_FixedPoint
 
    ! =================================================================
    ! Get functions
    ! =================================================================

     integer(c_int) function FSUNNonlinSolGetNumIters_FixedPoint(NLS, niters) &
         bind(C,name='SUNNonlinSolGetNumIters_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       integer(c_long)       :: niters
     end function FSUNNonlinSolGetNumIters_FixedPoint

     integer(c_int) function FSUNNonlinSolGetCurIter_FixedPoint(NLS, iter) &
         bind(C,name='SUNNonlinSolGetCurIter_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),   value :: NLS
       integer(c_int)       :: iter
     end function FSUNNonlinSolGetCurIter_FixedPoint
     
     integer(c_int) function FSUNNonlinSolGetSysFn_FixedPoint(NLS, SysFn) &
         bind(C,name='SUNNonlinSolGetSysFn_FixedPoint')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),   value :: NLS
       type(c_funptr)       :: SysFn
     end function FSUNNonlinSolGetSysFn_FixedPoint

   end interface

end module fsunnonlinsol_fixedpoint_mod
