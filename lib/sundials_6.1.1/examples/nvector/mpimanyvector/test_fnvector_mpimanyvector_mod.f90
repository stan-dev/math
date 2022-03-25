! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2022, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! -----------------------------------------------------------------
! This file tests the Fortran 2003 interface to the SUNDIALS
! manyvector N_Vector implementation.
! -----------------------------------------------------------------

module test_nvector_mpimanyvector
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_mpimanyvector_mod
  use fnvector_serial_mod
  use test_utilities
  implicit none
  include "mpif.h"

  integer(c_int), parameter  :: nsubvecs = 2
  integer(c_long), parameter :: N1       = 100        ! individual vector length
  integer(c_long), parameter :: N2       = 200        ! individual vector length
  integer(c_int),  parameter :: nv       = 3          ! length of vector arrays
  integer(c_long), parameter :: N        = N1 + N2    ! overall manyvector length
  integer(c_int), target     :: comm = MPI_COMM_WORLD ! default MPI communicator
  integer(c_int), pointer    :: commptr
  integer(c_int)             :: nprocs                ! number of MPI processes

contains

  integer function smoke_tests() result(ret)
    implicit none

    integer(c_long)         :: lenrw(1), leniw(1)     ! real and int work space size
    integer(c_long)         :: ival                   ! integer work value
    type(c_ptr)             :: cptr                   ! c_ptr work value
    real(c_double)          :: rval                   ! real work value
    real(c_double)          :: x1data(N1), x2data(N2) ! vector data array
    real(c_double), pointer :: xptr(:)                ! pointer to vector data array
    real(c_double)          :: nvarr(nv)              ! array of nv constants to go with vector array
    type(N_Vector), pointer :: x, y, z, tmp           ! N_Vectors
    type(c_ptr)             :: subvecs                ! MPIManyVector subvectors
    type(c_ptr)             :: xvecs, zvecs           ! C pointer to array of MPIManyVectors

    !===== Setup ====
    subvecs = FN_VNewVectorArray(nsubvecs)
    tmp  => FN_VMake_Serial(N1, x1data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 0, tmp)
    tmp  => FN_VMake_Serial(N2, x2data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 1, tmp)

    x => FN_VMake_MPIManyVector(MPI_COMM_WORLD, int(nsubvecs,8), subvecs, sunctx)
    call FN_VConst(ONE, x)
    y => FN_VClone_MPIManyVector(x)
    call FN_VConst(ONE, y)
    z => FN_VClone_MPIManyVector(x)
    call FN_VConst(ONE, z)

    xvecs = FN_VCloneVectorArray(nv, x)
    zvecs = FN_VCloneVectorArray(nv, z)
    nvarr = (/ ONE, ONE, ONE /)

    !===== Test =====

    ! test generic vector functions
    ival = FN_VGetVectorID_MPIManyVector(x)
    call FN_VSpace_MPIManyVector(x, lenrw, leniw)
    cptr = FN_VGetCommunicator(x)
    ival = FN_VGetLength_MPIManyVector(x)

    ! test standard vector operations
    call FN_VLinearSum_MPIManyVector(ONE, x, ONE, y, z)
    call FN_VConst_MPIManyVector(ONE, z)
    call FN_VProd_MPIManyVector(x, y, z)
    call FN_VDiv_MPIManyVector(x, y, z)
    call FN_VScale_MPIManyVector(ONE, x, y)
    call FN_VAbs_MPIManyVector(x, y)
    call FN_VInv_MPIManyVector(x, z)
    call FN_VAddConst_MPIManyVector(x, ONE, z)
    rval = FN_VDotProdLocal_MPIManyVector(x, y)
    rval = FN_VMaxNormLocal_MPIManyVector(x)
    rval = FN_VWrmsNorm_MPIManyVector(x, y)
    rval = FN_VWrmsNormMask_MPIManyVector(x, y, z)
    rval = FN_VMinLocal_MPIManyVector(x)
    rval = FN_VWL2Norm_MPIManyVector(x, y)
    rval = FN_VL1NormLocal_MPIManyVector(x)
    call FN_VCompare_MPIManyVector(ONE, x, y)
    ival = FN_VInvTestLocal_MPIManyVector(x, y)
    ival = FN_VConstrMaskLocal_MPIManyVector(z, x, y)
    rval = FN_VMinQuotientLocal_MPIManyVector(x, y)

    ! test fused vector operations
    ival = FN_VLinearCombination_MPIManyVector(nv, nvarr, xvecs, x)
    ival = FN_VScaleAddMulti_MPIManyVector(nv, nvarr, x, xvecs, zvecs)
    ival = FN_VDotProdMulti_MPIManyVector(nv, x, xvecs, nvarr)

    ! test vector array operations
    ival = FN_VLinearSumVectorArray_MPIManyVector(nv, ONE, xvecs, ONE, xvecs, zvecs)
    ival = FN_VScaleVectorArray_MPIManyVector(nv, nvarr, xvecs, zvecs)
    ival = FN_VConstVectorArray_MPIManyVector(nv, ONE, xvecs)
    ival = FN_VWrmsNormVectorArray_MPIManyVector(nv, xvecs, xvecs, nvarr)
    ival = FN_VWrmsNormMaskVectorArray_MPIManyVector(nv, xvecs, xvecs, x, nvarr)

    ! test the MPIManyVector specific operations
    ival = FN_VGetNumSubvectors_MPIManyVector(x)
    xptr => FN_VGetSubvectorArrayPointer_MPIManyVector(x, ival-1)
    ival = FN_VSetSubvectorArrayPointer_MPIManyVector(xptr, x, ival-1)
    tmp  => FN_VGetSubvector_MPIManyVector(x, ival-1)

    !==== Cleanup =====
    tmp => FN_VGetVecAtIndexVectorArray(subvecs, 0)
    call FN_VDestroy(tmp)
    tmp => FN_VGetVecAtIndexVectorArray(subvecs, 1)
    call FN_VDestroy(tmp)
    call FN_VDestroy_MPIManyVector(x)
    call FN_VDestroy_MPIManyVector(y)
    call FN_VDestroy_MPIManyVector(z)
    call FN_VDestroyVectorArray(xvecs, nv)
    call FN_VDestroyVectorArray(zvecs, nv)

    ret = 0

  end function smoke_tests

  integer function unit_tests() result(fails)
    use test_fnvector
    implicit none

    integer(c_int)           :: myid         ! the MPI rank
    real(c_double)           :: x1data(N1)   ! vector data array
    real(c_double)           :: x2data(N2)   ! vector data array
    type(N_Vector), pointer  :: x, tmp       ! N_Vectors
    type(c_ptr)              :: subvecs      ! subvectors of the MPIManyVector

    !===== Setup ====
    fails = 0

    call MPI_Comm_rank(comm, myid, fails)
    if (fails /= 0) then
      print *, '   FAILURE - MPI_COMM_RANK returned nonzero'
      stop 1
    endif

    subvecs = FN_VNewVectorArray(nsubvecs)
    tmp  => FN_VMake_Serial(N1, x1data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 0, tmp)
    tmp  => FN_VMake_Serial(N2, x2data, sunctx)
    call FN_VSetVecAtIndexVectorArray(subvecs, 1, tmp)

    x => FN_VMake_MPIManyVector(MPI_COMM_WORLD, int(nsubvecs,8), subvecs, sunctx)
    call FN_VConst(ONE, x)

    !==== tests ====
    fails = Test_FN_VMake(x, N, myid)
    fails = Test_FN_VGetArrayPointer(x, N, myid)
    fails = Test_FN_VLinearCombination(x, N, myid)

    !=== cleanup ====
    tmp => FN_VGetVecAtIndexVectorArray(subvecs, 0)
    call FN_VDestroy(tmp)
    tmp => FN_VGetVecAtIndexVectorArray(subvecs, 1)
    call FN_VDestroy(tmp)
    call FN_VDestroy_MPIManyVector(x)

  end function unit_tests

end module


integer(C_INT) function check_ans(ans, X, local_length) result(failure)
  use, intrinsic :: iso_c_binding
  use fnvector_mpimanyvector_mod
  use fsundials_nvector_mod
  use test_utilities
  implicit none

  real(C_DOUBLE)          :: ans
  type(N_Vector)          :: X
  type(N_Vector), pointer :: X0, X1
  integer(C_LONG)         :: local_length, i, x0len, x1len
  real(C_DOUBLE), pointer :: x0data(:), x1data(:)

  failure = 0

  X0 => FN_VGetSubvector_MPIManyVector(X, 0_8)
  X1 => FN_VGetSubvector_MPIManyVector(X, 1_8)
  x0len = FN_VGetLength(X0)
  x1len = FN_VGetLength(X1)
  x0data => FN_VGetArrayPointer(X0)
  x1data => FN_VGetArrayPointer(X1)

  if (local_length /= (x0len + x1len)) then
    failure = 1
    return
  endif

  do i = 1, x0len
    if (FNEQ(x0data(i), ans) > 0) then
      failure = failure + 1
    end if
  enddo

  do i = 1, x1len
    if (FNEQ(x1data(i), ans) > 0) then
      failure = failure + 1
    end if
  end do
end function check_ans


logical function has_data(X) result(failure)
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use test_utilities
  implicit none

  type(N_Vector)          :: X

  failure = .true.
end function has_data


program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_nvector_mpimanyvector

  !======== Declarations ========
  implicit none
  integer(c_int) :: fails = 0
  integer(c_int) :: myid

  !===== MPI Setup ====
  call MPI_Init(fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_INIT returned nonzero'
    stop 1
  endif

  call MPI_Comm_rank(comm, myid, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_COMM_RANK returned nonzero, proc', myid
    stop 1
  endif

  commptr => comm

  !============== Introduction =============
  if (myid == 0) print *, 'MPIManyVector N_Vector Fortran 2003 interface test'

  call Test_Init(c_loc(commptr))

  call MPI_Comm_size(comm, nprocs, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_COMM_SIZE returned nonzero, proc', myid
    stop 1
  endif

  fails = smoke_tests()
  if (fails /= 0) then
    print *, 'FAILURE: smoke tests failed, proc', myid
  else
    if (myid == 0) &
      print *, 'SUCCESS: all smoke tests passed'
  end if

  call MPI_BARRIER(comm, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_BARRIER returned nonzero, proc', myid
    stop 1
  endif

  fails = unit_tests()
  if (fails /= 0) then
    print *, 'FAILURE: n unit tests failed, proc', myid
    stop 1
  else
    if (myid == 0) print *,'    SUCCESS - all unit tests passed'
  end if

  call Test_Finalize()

  call MPI_Finalize(fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_FINALIZE returned nonzero, proc ', myid
    stop 1
  endif
end program main
