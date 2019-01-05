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
! the SUNDIALS Newton iteration nonlinear solver using the
! ISO_C_BINDING module.
! -----------------------------------------------------------------

module fsunnonlinsol_newton_mod

  !======= Interfaces =========
  interface

     ! =================================================================
     ! Constructors
     ! =================================================================

     type(c_ptr) function FSUNNonlinSol_Newton(y) &
         bind(C,name='SUNNonlinSol_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: y
     end function FSUNNonlinSol_Newton

     type(c_ptr) function FSUNNonlinSol_NewtonSens(cnt, y) &
         bind(C,name='SUNNonlinSol_NewtonSens')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: cnt
       type(c_ptr),    value :: y
     end function FSUNNonlinSol_NewtonSens

     ! =================================================================
     ! Destructors
     ! =================================================================

     integer(c_int) function FSUNNonlinSolFree_Newton(NLS) &
         bind(C,name='SUNNonlinSolFree_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolFree_Newton

     ! =================================================================
     ! Operations
     ! =================================================================

     integer(c_int) function FSUNNonlinSolGetType_Newton(NLS) &
         bind(C,name='SUNNonlinSolGetType_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolGetType_Newton

     integer(c_int) function FSUNNonlinSolInitialize_Newton(NLS) &
         bind(C,name='SUNNonlinSolInitialize_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: NLS
     end function FSUNNonlinSolInitialize_Newton

     integer(c_int) function FSUNNonlinSolSolve_Newton(NLS, y0, y, w, tol, &
                                                           callSetup, mem) &
         bind(C,name='SUNNonlinSolSolve_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_ptr),    value :: y0
       type(c_ptr),    value :: y
       type(c_ptr),    value :: w
       real(c_double), value :: tol
       integer(c_int), value :: callSetup
       type(c_ptr),    value :: mem
     end function FSUNNonlinSolSolve_Newton

    ! =================================================================
    ! Set functions
    ! =================================================================

     integer(c_int) function FSUNNonlinSolSetSysFn_Newton(NLS, SysFn) &
         bind(C,name='SUNNonlinSolSetSysFn_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: SysFn
     end function FSUNNonlinSolSetSysFn_Newton

     integer(c_int) function FSUNNonlinSolSetLSetupFn_Newton(NLS, LSetupFn) &
         bind(C,name='SUNNonlinSolSetLSetupFn_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: LSetupFn
     end function FSUNNonlinSolSetLSetupFn_Newton

     integer(c_int) function FSUNNonlinSolSetLSolveFn_Newton(NLS, LSolveFn) &
         bind(C,name='SUNNonlinSolSetLSolveFn_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: LSolveFn
     end function FSUNNonlinSolSetLSolveFn_Newton

     integer(c_int) function FSUNNonlinSolSetConvTestFn_Newton(NLS, CTestFN) &
         bind(C,name='SUNNonlinSolSetConvTestFn_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       type(c_funptr), value :: CTestFN
     end function FSUNNonlinSolSetConvTestFn_Newton

     integer(c_int) function FSUNNonlinSolSetMaxIters_Newton(NLS, maxiters) &
         bind(C,name='SUNNonlinSolSetMaxIters_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       integer(c_int), value :: maxiters
     end function FSUNNonlinSolSetMaxIters_Newton

    ! =================================================================
    ! Get functions
    ! =================================================================

     integer(c_int) function FSUNNonlinSolGetNumIters_Newton(NLS, niters) &
         bind(C,name='SUNNonlinSolGetNumIters_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: NLS
       integer(c_long)       :: niters
     end function FSUNNonlinSolGetNumIters_Newton

     integer(c_int) function FSUNNonlinSolGetCurIter_Newton(NLS, iter) &
         bind(C,name='SUNNonlinSolGetCurIter_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),   value :: NLS
       integer(c_int)       :: iter
     end function FSUNNonlinSolGetCurIter_Newton

     integer(c_int) function FSUNNonlinSolGetSysFn_Newton(NLS, SysFn) &
         bind(C,name='SUNNonlinSolGetSysFn_Newton')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),   value :: NLS
       type(c_funptr)       :: SysFn
     end function FSUNNonlinSolGetSysFn_Newton

   end interface

end module fsunnonlinsol_newton_mod
