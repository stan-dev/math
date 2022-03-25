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
! parallel N_Vector implementation.
! -----------------------------------------------------------------

module test_nvector_parallel
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use fnvector_parallel_mod
  use test_utilities
  implicit none
  include "mpif.h"

  integer(c_long), parameter :: local_length = 100    ! vector local length
  integer(c_int),  parameter :: nv = 3                ! length of vector arrays
  integer(c_int),  parameter :: ns = 2                ! number of vector arrays

  integer(c_int), target  :: comm = MPI_COMM_WORLD ! default MPI communicator
  integer(c_int), pointer :: commptr
  integer(c_long)         :: global_length ! vector global_length
  integer(c_int)          :: nprocs        ! number of MPI processes
  contains

  integer function smoke_tests() result(ret)
    implicit none

    integer(c_long)         :: lenrw(1), leniw(1)  ! real and int work space size
    integer(c_long)         :: ival                ! integer work value
    type(c_ptr)             :: cptr                ! c_ptr work value
    real(c_double)          :: rval                ! real work value
    real(c_double)          :: xdata(local_length) ! vector data array
    real(c_double), pointer :: xptr(:)             ! pointer to vector data array
    real(c_double)          :: nvarr(nv)           ! array of nv constants to go with vector array
    type(N_Vector), pointer :: x, y, z, tmp        ! N_Vectors
    type(c_ptr)             :: xvecs, zvecs        ! C pointer to array of C pointers to N_Vectors

    !===== Setup ====
    x => FN_VMake_Parallel(comm, local_length, global_length, xdata, sunctx)
    call FN_VConst(ONE, x)
    y => FN_VClone_Parallel(x)
    call FN_VConst(ONE, y)
    z => FN_VClone_Parallel(x)
    call FN_VConst(ONE, z)

    xvecs = FN_VCloneVectorArray(nv, x)
    zvecs = FN_VCloneVectorArray(nv, z)
    nvarr = (/ ONE, ONE, ONE /)

    !===== Test =====

    ! test constructors
    tmp => FN_VNewEmpty_Parallel(comm, local_length, global_length, sunctx)
    call FN_VDestroy_Parallel(tmp)
    tmp => FN_VMake_Parallel(comm, local_length, global_length, xdata, sunctx)
    call FN_VDestroy_Parallel(tmp)
    tmp => FN_VNew_Parallel(comm, local_length, global_length, sunctx)
    call FN_VDestroy_Parallel(tmp)
    tmp => FN_VCloneEmpty_Parallel(x)
    call FN_VDestroy_Parallel(tmp)

    ! test generic vector functions
    ival = FN_VGetVectorID_Parallel(x)
    call FN_VSpace_Parallel(x, lenrw, leniw)
    xptr => FN_VGetArrayPointer_Parallel(x)
    call FN_VSetArrayPointer_Parallel(xdata, x)
    cptr = FN_VGetCommunicator_Parallel(x)
    ival = FN_VGetLength_Parallel(x)

    ! test standard vector operations
    call FN_VLinearSum_Parallel(ONE, x, ONE, y, z)
    call FN_VConst_Parallel(ONE, z)
    call FN_VProd_Parallel(x, y, z)
    call FN_VDiv_Parallel(x, y, z)
    call FN_VScale_Parallel(ONE, x, y)
    call FN_VAbs_Parallel(x, y)
    call FN_VInv_Parallel(x, z)
    call FN_VAddConst_Parallel(x, ONE, z)
    rval = FN_VDotProd_Parallel(x, y)
    rval = FN_VMaxNorm_Parallel(x)
    rval = FN_VWrmsNorm_Parallel(x, y)
    rval = FN_VWrmsNormMask_Parallel(x, y, z)
    rval = FN_VMin_Parallel(x)
    rval = FN_VWL2Norm_Parallel(x, y)
    rval = FN_VL1Norm_Parallel(x)
    call FN_VCompare_Parallel(ONE, x, y)
    ival = FN_VInvTest_Parallel(x, y)
    ival = FN_VConstrMask_Parallel(z, x, y)
    rval = FN_VMinQuotient_Parallel(x, y)

    ! test fused vector operations
    ival = FN_VLinearCombination_Parallel(nv, nvarr, xvecs, x)
    ival = FN_VScaleAddMulti_Parallel(nv, nvarr, x, xvecs, zvecs)
    ival = FN_VDotProdMulti_Parallel(nv, x, xvecs, nvarr)

    ! test vector array operations
    ival = FN_VLinearSumVectorArray_Parallel(nv, ONE, xvecs, ONE, xvecs, zvecs)
    ival = FN_VScaleVectorArray_Parallel(nv, nvarr, xvecs, zvecs)
    ival = FN_VConstVectorArray_Parallel(nv, ONE, xvecs)
    ival = FN_VWrmsNormVectorArray_Parallel(nv, xvecs, xvecs, nvarr)
    ival = FN_VWrmsNormMaskVectorArray_Parallel(nv, xvecs, xvecs, x, nvarr)

    !==== Cleanup =====
    call FN_VDestroy_Parallel(x)
    call FN_VDestroy_Parallel(y)
    call FN_VDestroy_Parallel(z)
    call FN_VDestroyVectorArray(xvecs, nv)
    call FN_VDestroyVectorArray(zvecs, nv)

    ret = 0

  end function smoke_tests

  integer function unit_tests() result(fails)
    use test_fnvector
    implicit none

    ! local variables
    integer(c_int)  :: myid                ! the MPI rank
    real(c_double)  :: xdata(local_length) ! vector data array
    type(N_Vector), pointer  :: x          ! N_Vectors

    !===== Setup ====
    fails = 0

    call MPI_Comm_rank(comm, myid, fails)
    if (fails /= 0) then
      print *, '   FAILURE - MPI_COMM_RANK returned nonzero'
      stop 1
    endif

    x => FN_VMake_Parallel(comm, local_length, global_length, xdata, sunctx)
    call FN_VConst(ONE, x)

    !==== tests ====
    fails = Test_FN_VMake(x, local_length, myid)
    fails = Test_FN_VGetArrayPointer(x, local_length, myid)
    fails = Test_FN_VLinearCombination(x, local_length, myid)

    !=== cleanup ====
    call FN_VDestroy_Parallel(x)

  end function unit_tests

end module


integer(C_INT) function check_ans(ans, X, local_length) result(failure)
  use, intrinsic :: iso_c_binding
  use fsundials_nvector_mod
  use test_utilities
  implicit none

  real(C_DOUBLE)          :: ans
  type(N_Vector)          :: X
  integer(C_LONG)         :: local_length, i
  real(C_DOUBLE), pointer :: Xdata(:)

  failure = 0

  Xdata => FN_VGetArrayPointer(X)
  do i = 1, local_length
    if (FNEQ(Xdata(i), ans) > 0) then
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
  real(C_DOUBLE), pointer :: xptr(:)

  xptr => FN_VGetArrayPointer(x)
  failure = associated(xptr)
end function has_data


program main
  !======== Inclusions ==========
  use, intrinsic :: iso_c_binding
  use test_nvector_parallel

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
  if (myid == 0) print *, 'Parallel N_Vector Fortran 2003 interface test'

  call Test_Init(c_loc(commptr))

  call MPI_Comm_size(comm, nprocs, fails)
  if (fails /= 0) then
    print *, 'FAILURE: MPI_COMM_SIZE returned nonzero, proc', myid
    stop 1
  endif
  global_length = nprocs*local_length

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
