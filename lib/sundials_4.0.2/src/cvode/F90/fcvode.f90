! ------------------------------------------------------------------
! Programmer(s): David J. Gardner @ LLNL
!                Daniel R. Reynolds @ SMU
! ------------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! This file contains a Fortran module for interfacing directly with
! CVODE using the ISO_C_BINDING module.
! ------------------------------------------------------------------

module fcvode_mod

  use, intrinsic :: iso_c_binding, only : c_int

  ! =================================================================
  !              C V O D E     C O N S T A N T S
  ! =================================================================

  ! -----------------------------------------------------------------
  ! Enumerations for inputs to CVodeCreate and CVode.
  ! -----------------------------------------------------------------

  ! lmm
  integer(c_int), parameter :: CV_ADAMS = 1
  integer(c_int), parameter :: CV_BDF   = 2

  ! itask
  integer(c_int), parameter :: CV_NORMAL   = 1
  integer(c_int), parameter :: CV_ONE_STEP = 2

  ! -----------------------------------------------------------------
  ! CVODE return flags
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CV_SUCCESS           =   0
  integer(c_int), parameter :: CV_TSTOP_RETURN      =   1
  integer(c_int), parameter :: CV_ROOT_RETURN       =   2

  integer(c_int), parameter :: CV_WARNING           =  99

  integer(c_int), parameter :: CV_TOO_MUCH_WORK     =  -1
  integer(c_int), parameter :: CV_TOO_MUCH_ACC      =  -2
  integer(c_int), parameter :: CV_ERR_FAILURE       =  -3
  integer(c_int), parameter :: CV_CONV_FAILURE      =  -4

  integer(c_int), parameter :: CV_LINIT_FAIL        =  -5
  integer(c_int), parameter :: CV_LSETUP_FAIL       =  -6
  integer(c_int), parameter :: CV_LSOLVE_FAIL       =  -7
  integer(c_int), parameter :: CV_RHSFUNC_FAIL      =  -8
  integer(c_int), parameter :: CV_FIRST_RHSFUNC_ERR =  -9
  integer(c_int), parameter :: CV_REPTD_RHSFUNC_ERR =  -10
  integer(c_int), parameter :: CV_UNREC_RHSFUNC_ERR =  -11
  integer(c_int), parameter :: CV_RTFUNC_FAIL       =  -12
  integer(c_int), parameter :: CV_NLS_INIT_FAIL     =  -13
  integer(c_int), parameter :: CV_NLS_SETUP_FAIL    =  -14
  integer(c_int), parameter :: CV_CONSTR_FAIL       =  -15

  integer(c_int), parameter :: CV_MEM_FAIL          =  -20
  integer(c_int), parameter :: CV_MEM_NULL          =  -21
  integer(c_int), parameter :: CV_ILL_INPUT         =  -22
  integer(c_int), parameter :: CV_NO_MALLOC         =  -23
  integer(c_int), parameter :: CV_BAD_K             =  -24
  integer(c_int), parameter :: CV_BAD_T             =  -25
  integer(c_int), parameter :: CV_BAD_DKY           =  -26
  integer(c_int), parameter :: CV_TOO_CLOSE         =  -27
  integer(c_int), parameter :: CV_VECTOROP_ERR      =  -28

  ! -----------------------------------------------------------------
  ! CVLS return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CVLS_SUCCESS         =  0
  integer(c_int), parameter :: CVLS_MEM_NULL        = -1
  integer(c_int), parameter :: CVLS_LMEM_NULL       = -2
  integer(c_int), parameter :: CVLS_ILL_INPUT       = -3
  integer(c_int), parameter :: CVLS_MEM_FAIL        = -4
  integer(c_int), parameter :: CVLS_PMEM_NULL       = -5
  integer(c_int), parameter :: CVLS_JACFUN_UNRECVR  = -6
  integer(c_int), parameter :: CVLS_JACFUNC_RECVR   = -7
  integer(c_int), parameter :: CVLS_SUNMAT_FAIL     = -8
  integer(c_int), parameter :: CVLS_SUNLS_FAIL      = -9

  ! -----------------------------------------------------------------
  ! CVDIAG return values
  ! -----------------------------------------------------------------

  integer(c_int), parameter :: CVDIAG_SUCCESS         =  0
  integer(c_int), parameter :: CVDIAG_MEM_NULL        = -1
  integer(c_int), parameter :: CVDIAG_LMEM_NULL       = -2
  integer(c_int), parameter :: CVDIAG_ILL_INPUT       = -3
  integer(c_int), parameter :: CVDIAG_MEM_FAIL        = -4

  ! Additional last_flag values
  integer(c_int), parameter :: CVDIAG_INV_FAIL        = -5
  integer(c_int), parameter :: CVDIAG_RHSFUNC_UNRECVR = -6
  integer(c_int), parameter :: CVDIAG_RHSFUNC_RECVR   = -7

  ! =================================================================
  !          U S E R - C A L L A B L E   R O U T I N E S
  ! =================================================================

  interface

     ! =================================================================
     ! Interfaces for cvode.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! CVodeCreate
     ! -----------------------------------------------------------------

     type(c_ptr) function FCVodeCreate(lmm) &
          bind(C,name='CVodeCreate')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: lmm
     end function FCVodeCreate

     ! -----------------------------------------------------------------
     ! Integrator optional input specification functions
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeSetErrHandlerFn(cvode_mem, ehfun, eh_data) &
          bind(C,name='CVodeSetErrHandlerFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: ehfun
       type(c_ptr),    value :: eh_data
     end function FCVodeSetErrHandlerFn

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: CVodeSetErrFile
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeSetUserData(cvode_mem, user_data) &
          bind(C,name='CVodeSetUserData')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: user_data
     end function FCVodeSetUserData

     integer(c_int) function FCVodeSetMaxOrd(cvode_mem, maxord) &
          bind(C,name='CVodeSetMaxOrd')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxord
     end function FCVodeSetMaxOrd

     integer(c_int) function FCVodeSetMaxNumSteps(cvode_mem, mxsteps) &
          bind(C,name='CVodeSetMaxNumSteps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: mxsteps
     end function FCVodeSetMaxNumSteps

     integer(c_int) function FCVodeSetMaxHnilWarns(cvode_mem, mxhnil) &
          bind(C,name='CVodeSetMaxHnilWarns')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: mxhnil
     end function FCVodeSetMaxHnilWarns

     integer(c_int) function FCVodeSetStabLimDet(cvode_mem, stldet) &
          bind(C,name='CVodeSetStabLimDet')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: stldet
     end function FCVodeSetStabLimDet

     integer(c_int) function FCVodeSetInitStep(cvode_mem, hin) &
          bind(C,name='CVodeSetInitStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: hin
     end function FCVodeSetInitStep

     integer(c_int) function FCVodeSetMinStep(cvode_mem, hmin) &
          bind(C,name='CVodeSetMinStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: hmin
     end function FCVodeSetMinStep

     integer(c_int) function FCVodeSetMaxStep(cvode_mem, hmax) &
          bind(C,name='CVodeSetMaxStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: hmax
     end function FCVodeSetMaxStep

     integer(c_int) function FCVodeSetStopTime(cvode_mem, tstop) &
          bind(C,name='CVodeSetStopTime')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: tstop
     end function FCVodeSetStopTime

     integer(c_int) function FCVodeSetMaxErrTestFails(cvode_mem, maxnef) &
          bind(C,name='CVodeSetMaxErrTestFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxnef
     end function FCVodeSetMaxErrTestFails

     integer(c_int) function FCVodeSetMaxNonlinIters(cvode_mem, maxcor) &
          bind(C,name='CVodeSetMaxNonlinIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxcor
     end function FCVodeSetMaxNonlinIters

     integer(c_int) function FCVodeSetMaxConvFails(cvode_mem, maxncf) &
          bind(C,name='CVodeSetMaxConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: maxncf
     end function FCVodeSetMaxConvFails

     integer(c_int) function FCVodeSetNonlinConvCoef(cvode_mem, nlscoef) &
          bind(C,name='CVodeSetNonlinConvCoef')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: nlscoef
     end function FCVodeSetNonlinConvCoef

     integer(c_int) function FCVodeSetConstraints(cvode_mem, constraints) &
         bind(C,name='CVodeSetConstraints')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: constraints
     end function FCVodeSetConstraints

     integer(c_int) function FCVodeSetIterType(cvode_mem, iter) &
          bind(C,name='CVodeSetIterType')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: iter
     end function FCVodeSetIterType

     integer(c_int) function FCVodeSetRootDirection(cvode_mem, rootdir) &
          bind(C,name='CVodeSetRootDirection')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: rootdir
     end function FCVodeSetRootDirection

     integer(c_int) function FCVodeSetNoInactiveRootWarn(cvode_mem) &
         bind(C,name='CVodeSetNoInactiveRootWarn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
     end function FCVodeSetNoInactiveRootWarn

     integer(c_int) function FCVodeSetNonlinearSolver(cvode_mem, NLS) &
        bind(C,name='CVodeSetNonlinearSolver')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: NLS
     end function FCVodeSetNonlinearSolver

     ! -----------------------------------------------------------------
     ! CVodeInit
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeInit(cvode_mem, f, t0, y0) &
          bind(C,name='CVodeInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: f
       real(c_double), value :: t0
       type(c_ptr),    value :: y0
     end function FCVodeInit

     ! -----------------------------------------------------------------
     ! CVodeReInit
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeReInit(cvode_mem, t0, y0) &
          bind(C,name='CVodeReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: t0
       type(c_ptr),    value :: y0
     end function FCVodeReInit

     ! -----------------------------------------------------------------
     ! CVodeSStolerances
     ! CVodeSVtolerances
     ! CVodeWFtolerances
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeSStolerances(cvode_mem, reltol, abstol) &
          bind(C,name='CVodeSStolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: reltol
       real(c_double), value :: abstol
     end function FCVodeSStolerances

     integer(c_int) function FCVodeSVtolerances(cvode_mem, reltol, abstol) &
          bind(C,name='CVodeSVtolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: reltol
       type(c_ptr),    value :: abstol
     end function FCVodeSVtolerances

     integer(c_int) function FCVodeWFtolerances(cvode_mem, efun) &
          bind(C,name='CVodeWFtolerances')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: efun
     end function FCVodeWFtolerances

     ! -----------------------------------------------------------------
     ! CVodeRootInit
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeRootInit(cvode_mem, nrtfn, g) &
          bind(C,name='CVodeRootInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       integer(c_int), value :: nrtfn
       type(c_funptr), value :: g
     end function FCVodeRootInit

     ! -----------------------------------------------------------------
     ! CVode
     ! -----------------------------------------------------------------

     integer(c_int) function FCVode(cvode_mem, tout, yout, tret, itask) &
          bind(C,name='CVode')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: tout
       type(c_ptr),    value :: yout
       real(c_double)        :: tret
       integer(c_int), value :: itask
     end function FCVode

     ! -----------------------------------------------------------------
     ! CVodeGetDky
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetDky(cvode_mem, t, k, dky) &
          bind(C,name='CVodeGetDky')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: t
       integer(c_int), value :: k
       type(c_ptr),    value :: dky
     end function FCVodeGetDky

     ! -----------------------------------------------------------------
     ! Integrator optional output extraction functions
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetWorkSpace(cvode_mem, lenrw, leniw) &
          bind(C,name='CVodeGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrw
       integer(c_long)    :: leniw
     end function FCVodeGetWorkSpace

     integer(c_int) function FCVodeGetNumSteps(cvode_mem, nsteps) &
          bind(C,name='CVodeGetNumSteps')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nsteps
     end function FCVodeGetNumSteps

     integer(c_int) function FCVodeGetNumRhsEvals(cvode_mem, nfevals) &
          bind(C,name='CVodeGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevals
     end function FCVodeGetNumRhsEvals

     integer(c_int) function FCVodeGetNumLinSolvSetups(cvode_mem, nlinsetups) &
          bind(C,name='CVodeGetNumLinSolvSetups')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nlinsetups
     end function FCVodeGetNumLinSolvSetups

     integer(c_int) function FCVodeGetNumErrTestFails(cvode_mem, netfails) &
          bind(C,name='CVodeGetNumErrTestFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: netfails
     end function FCVodeGetNumErrTestFails

     integer(c_int) function FCVodeGetLastOrder(cvode_mem, qlast) &
          bind(C,name='CVodeGetLastOrder')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_int)    :: qlast
     end function FCVodeGetLastOrder

     integer(c_int) function FCVodeGetCurrentOrder(cvode_mem, qcur) &
          bind(C,name='CVodeGetCurrentOrder')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_int)     :: qcur
     end function FCVodeGetCurrentOrder

     integer(c_int) function FCVodeGetNumStabLimOrderReds(cvode_mem, nslred) &
          bind(C,name='CVodeGetNumStabLimOrderReds')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nslred
     end function FCVodeGetNumStabLimOrderReds

     integer(c_int) function FCVodeGetActualInitStep(cvode_mem, hinused) &
          bind(C,name='CVodeGetActualInitStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: hinused
     end function FCVodeGetActualInitStep

     integer(c_int) function FCVodeGetLastStep(cvode_mem, hlast) &
          bind(C,name='CVodeGetLastStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: hlast
     end function FCVodeGetLastStep

     integer(c_int) function FCVodeGetCurrentStep(cvode_mem, hcur) &
         bind(C,name='CVodeGetCurrentStep')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: hcur
     end function FCVodeGetCurrentStep

     integer(c_int) function FCVodeGetCurrentTime(cvode_mem, tcur) &
          bind(C,name='CVodeGetCurrentTime')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: tcur
     end function FCVodeGetCurrentTime

     integer(c_int) function FCVodeGetTolScaleFactor(cvode_mem, tolsfac) &
          bind(C,name='CVodeGetTolScaleFactor')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       real(c_double)     :: tolsfac
     end function FCVodeGetTolScaleFactor

     integer(c_int) function FCVodeGetErrWeights(cvode_mem, eweight) &
          bind(C,name='CVodeGetEffWeights')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: eweight
     end function FCVodeGetErrWeights

     integer(c_int) function FCVodeGetEstLocalErrors(cvode_mem, ele) &
          bind(C,name='CVodeGetEstLocalErrors')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: ele
     end function FCVodeGetEstLocalErrors

     integer(c_int) function FCVodeGetNumGEvals(cvode_mem, ngevals) &
          bind(C,name='CVodeGetNumGEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: ngevals
     end function FCVodeGetNumGEvals

     integer(c_int) function FCVodeGetRootInfo(cvode_mem, rootsfound) &
          bind(C,name='CVodeGetRootInfo')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_int)     :: rootsfound
     end function FCVodeGetRootInfo

     integer(c_int) function FCVodeGetIntegratorStats(cvode_mem, nsteps, nfevals, &
          nlinsetups, netfails, qlast, qcur, hinused, hlast, hcur, tcur) &
          bind(C,name='CVodeGetIntegratorStats')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nsteps
       integer(c_long)    :: nfevals
       integer(c_long)    :: nlinsetups
       integer(c_long)    :: netfails
       integer(c_int)     :: qlast
       integer(c_int)     :: qcur
       real(c_double)     :: hinused
       real(c_double)     :: hlast
       real(c_double)     :: hcur
       real(c_double)     :: tcur
     end function FCVodeGetIntegratorStats

     ! -----------------------------------------------------------------
     ! Nonlinear solver optional output extraction functions
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetNumNonlinSolvIters(cvode_mem, nniters) &
          bind(C,name='CVodeGetNumNonlinSolvIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nniters
     end function FCVodeGetNumNonlinSolvIters

     integer(c_int) function FCVodeGetNumNonlinSolvConvFails(cvode_mem, nncfails) &
          bind(C,name='CVodeGetNumNonlinSolvConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nncfails
     end function FCVodeGetNumNonlinSolvConvFails

     integer(c_int) function FCVodeGetNonlinSolvStats(cvode_mem, nniters, nncfails) &
          bind(C,name='CVodeGetNonlinSolvStats')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nniters
       integer(c_long)    :: nncfails
     end function FCVodeGetNonlinSolvStats

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: CVodeGetReturnFlagName
     ! -----------------------------------------------------------------

     ! -----------------------------------------------------------------
     ! CVodeFree
     ! -----------------------------------------------------------------

     subroutine FCVodeFree(cvode_mem) &
          bind(C,name='CVodeFree')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: cvode_mem ! DO NOT use value attribute input is void**
     end subroutine FCVodeFree

     ! =================================================================
     ! Interfaces from cvode_ls.h
     ! =================================================================

     integer(c_int) function FCVodeSetLinearSolver(cvode_mem, LS, A) &
          bind(C,name='CVodeSetLinearSolver')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       type(c_ptr), value :: LS
       type(c_ptr), value :: A
     end function FCVodeSetLinearSolver

     ! -----------------------------------------------------------------
     ! Optional inputs to the CVLS linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeSetJacFn(cvode_mem, jac) &
          bind(C,name='CVodeSetJacFn')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: jac
     end function FCVodeSetJacFn

     integer(c_int) function FCVodeSetMaxStepsBetweenJac(cvode_mem, msbj) &
          bind(C,name='CVodeSetMaxStepsBetweenJac')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: msbj
     end function FCVodeSetMaxStepsBetweenJac
     
     integer(c_int) function FCVodeSetEpsLin(cvode_mem, eplifac) &
          bind(C,name='CVodeSetEpsLin')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       real(c_double), value :: eplifac
     end function FCVodeSetEpsLin

     integer(c_int) function FCVodeSetPreconditioner(cvode_mem, pset, psolve) &
          bind(C,name='CVodeSetPreconditioner')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: pset
       type(c_funptr), value :: psolve
     end function FCVodeSetPreconditioner

     integer(c_int) function FCVodeSetJacTimes(cvode_mem, jtsetup, jtimes) &
          bind(C,name='CVodeSetJacTimes')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),    value :: cvode_mem
       type(c_funptr), value :: jtsetup
       type(c_funptr), value :: jtimes
     end function FCVodeSetJacTimes

     ! -----------------------------------------------------------------
     ! Optional outputs from the CVLS linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FCVodeGetLinWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVodeGetLinWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVodeGetLinWorkSpace
     
     integer(c_int) function FCVodeGetNumJacEvals(cvode_mem, njevals) &
          bind(C,name='CVodeGetNumJacEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: njevals
     end function FCVodeGetNumJacEvals

     integer(c_int) function FCVodeGetNumPrecEvals(cvode_mem, npevals) &
         bind(C,name='CVodeGetNumPrecEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: npevals
     end function FCVodeGetNumPrecEvals

     integer(c_int) function FCVodeGetNumPrecSolves(cvode_mem, npsolves) &
         bind(C,name='CVodeGetNumPrecSolves')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: npsolves
     end function FCVodeGetNumPrecSolves

     integer(c_int) function FCVodeGetNumLinIters(cvode_mem, nliters) &
         bind(C,name='CVodeGetNumLinIters')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nliters
     end function FCVodeGetNumLinIters

     integer(c_int) function FCVodeGetNumLinConvFails(cvode_mem, nlcfails) &
         bind(C,name='CVodeGetNumLinConvFails')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nlcfails
     end function FCVodeGetNumLinConvFails

     integer(c_int) function FCVodeGetNumJTSetupEvals(cvode_mem, njtsetups) &
         bind(C,name='CVodeGetNumJTSetupEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: njtsetups
     end function FCVodeGetNumJTSetupEvals

     integer(c_int) function FCVodeGetNumJtimesEvals(cvode_mem, njvevals) &
         bind(C,name='CVodeGetNumJtimesEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: njvevals
     end function FCVodeGetNumJtimesEvals

     integer(c_int) function FCVodeGetNumLinRhsEvals(cvode_mem, nfevalsLS) &
          bind(C,name='CVodeGetNumLinRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsLS
     end function FCVodeGetNumLinRhsEvals

     integer(c_int) function FCVodeGetLastLinFlag(cvode_mem, flag) &
          bind(C,name='CVodeGetLastLinFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: flag
     end function FCVodeGetLastLinFlag

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: CVodeGetLinReturnFlagName
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Interfaces for cvode_diag.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! CVDiag
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDiag(cvode_mem) &
          bind(C,name='CVDiag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
     end function FCVDiag

     ! -----------------------------------------------------------------
     ! Optional outputs from the CVDIAG linear solver
     ! -----------------------------------------------------------------

     integer(c_int) function FCVDiagGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVDiagGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVDiagGetWorkSpace

     integer(c_int) function FCVDiagGetNumRhsEvals(cvode_mem, nfevalsLS) &
          bind(C,name='CVDiagGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsLS
     end function FCVDiagGetNumRhsEvals

     integer(c_int) function FCVDiagGetLastFlag(cvode_mem, flag) &
          bind(C,name='CVDiagGetLastFlag')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: flag
     end function FCVDiagGetLastFlag

     ! -----------------------------------------------------------------
     ! NOT INTERFACED: CVDiagGetReturnFlagName
     ! -----------------------------------------------------------------

     ! =================================================================
     ! Interfaces from cvode_bandpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! CVBandPrecInit
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBandPrecInit(cvode_mem, N, mu, ml) &
          bind(C,name='CVBandPrecInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: N
       integer(c_long), value :: mu
       integer(c_long), value :: ml
     end function FCVBandPrecInit

     ! -----------------------------------------------------------------
     ! Optional output functions : CVBandPrecGet*
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBandPrecGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVBandPrecGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVBandPrecGetWorkSpace

     integer(c_int) function FCVBandPrecGetNumRhsEvals(cvode_mem, nfevalsBP) &
          bind(C,name='CVBandPrecGetNumRhsEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: nfevalsBP
     end function FCVBandPrecGetNumRhsEvals

     ! =================================================================
     ! Interfaces for cvode_bbdpre.h
     ! =================================================================

     ! -----------------------------------------------------------------
     ! CVBBDPrecInit
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBBDPrecInit(cvode_mem, Nlocal, mudq, mldq, &
          mukeep, mlkeep, dqrely, gloc, cfn) &
          bind(C,name='CVBBDPrecInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: Nlocal
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       integer(c_long), value :: mukeep
       integer(c_long), value :: mlkeep
       real(c_double),  value :: dqrely
       type(c_funptr),  value :: gloc
       type(c_funptr),  value :: cfn
     end function FCVBBDPrecInit

     ! -----------------------------------------------------------------
     ! CVBBDPrecReInit
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBBDPrecReInit(cvode_mem, mudq, mldq, dqrely) &
          bind(C,name='CVBBNPrecReInit')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr),     value :: cvode_mem
       integer(c_long), value :: mudq
       integer(c_long), value :: mldq
       real(c_double),  value :: dqrely
     end function FCVBBDPrecReInit

     ! -----------------------------------------------------------------
     ! BBDPRE optional output extraction routines
     ! -----------------------------------------------------------------

     integer(c_int) function FCVBBDPrecGetWorkSpace(cvode_mem, lenrwLS, leniwLS) &
          bind(C,name='CVBBDPrecGetWorkSpace')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: lenrwLS
       integer(c_long)    :: leniwLS
     end function FCVBBDPrecGetWorkSpace

     integer(c_int) function FCVBBDPrecGetNumGfnEvals(cvode_mem, ngevalsBBDP) &
          bind(C,name='CVBBDPrecGetNumGfnEvals')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: cvode_mem
       integer(c_long)    :: ngevalsBBDP
     end function FCVBBDPrecGetNumGfnEvals

  end interface

end module fcvode_mod
