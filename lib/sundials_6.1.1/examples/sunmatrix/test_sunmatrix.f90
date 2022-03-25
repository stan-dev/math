! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! Acknowledgements: These testing routines are based on
! test_sunmatrix.c written by David Gardner @ LLNL and Daniel
! R. Reynolds @ SMU.
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
! These test functions are designed to check the SWIG generated
! Fortran interface to a SUNMatrix module implementation.
! -----------------------------------------------------------------

module test_sunmatrix
  use, intrinsic :: iso_c_binding
  use test_utilities
  use fsundials_nvector_mod
  use fsundials_matrix_mod
  use fsundials_types_mod

  implicit none

  logical, parameter :: print_all_ranks = .false.

  ! functions implemented in specific matrix tests
  integer(C_INT), external :: check_matrix
  integer(C_INT), external :: check_matrix_entry
  logical,        external :: is_square

contains

  subroutine TEST_STATUS(frmt, myrank)
    use, intrinsic :: iso_c_binding

    implicit none

    character(LEN=*) :: frmt
    integer(C_INT)   :: myrank

    if (print_all_ranks) then
      write(*,'(A,I0,A,A)') 'process ', myrank, ': ', frmt
    else
      write(*,*) frmt
    end if

  end subroutine TEST_STATUS

  subroutine TEST_STATUS2(frmt, retval, myrank)
    use, intrinsic :: iso_c_binding

    implicit none

    character(LEN=*) :: frmt
    integer(C_INT)   :: myrank
    integer(C_INT)   :: retval

    if (print_all_ranks) then
      write(*,'(A,I0,A,A,I0)') 'process ', myrank, ': ', frmt, retval
    else
      write(*,'(A,I0)') frmt, retval
    end if

  end subroutine TEST_STATUS2

  integer(C_INT) function check_vector(x, y, tol) result(failure)
    use, intrinsic :: iso_c_binding

    implicit none

    type(N_Vector)  :: x, y
    real(C_DOUBLE)  :: tol
    integer(C_LONG) :: i, xlen, ylen
    real(C_DOUBLE), pointer :: xdata(:), ydata(:)

    failure = 0

    xdata => FN_VGetArrayPointer(x)
    ydata => FN_VGetArrayPointer(y)

    xlen = FN_VGetLength(x)
    ylen = FN_VGetLength(y)

    if (xlen /= ylen) then
      print *, 'ERROR: check_vector: different data array lengths'
      failure = 1
      return
    end if

    do i = 1, xlen
      failure = failure + FNEQTOL(xdata(i), ydata(i), tol)
    end do

  end function check_vector

  integer(C_INT) function Test_FSUNMatGetID(A, sunid, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    type(SUNMatrix)       :: A
    integer(SUNMatrix_ID) :: sunid, mysunid
    integer(C_INT)        :: myid

    failure = 0

    mysunid = FSUNMatGetID(A)

    if (sunid /= mysunid) then
      call TEST_STATUS(">>> FAILED test -- SUNMatGetID ", myid)
      failure = 1
    else
      call TEST_STATUS("    PASSED test -- SUNMatGetID ", myid)
    end if

  end function Test_FSUNMatGetID

  ! --------------------------------------------------------------------
  ! SUNMatClone Test
  ! NOTE: This routine depends on SUNMatCopy to check matrix data.
  ! --------------------------------------------------------------------
  integer(C_INT) function Test_FSUNMatClone(A, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    type(SUNMatrix) :: A
    real(C_DOUBLE)  :: tol = 10*UNIT_ROUNDOFF
    type(SUNMatrix), pointer ::  B

    failure = 0

    ! clone vector
    B => FSUNMatClone(A)

    ! check cloned matrix
    if (.not. associated(B)) then
      call TEST_STATUS(">>> FAILED test -- SUNMatClone (after SUNMatClone, B == NULL)", myid)
      failure = 1
      return
    end if

    failure = FSUNMatCopy(A, B)
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- FSUNMatClone (error using SUNMatCopy) ", myid)
      call FSUNMatDestroy(B)
      return
    end if

    failure = check_matrix(B, A, tol)
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatClone (failure at check_matrix)", myid)
      call FSUNMatDestroy(B)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatClone ", myid)
    end if

    call FSUNMatDestroy(B)

  end function Test_FSUNMatClone

  integer(C_INT) function Test_FSUNMatZero(A, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    real(C_DOUBLE)  :: tol = 10*UNIT_ROUNDOFF
    type(SUNMatrix) :: A
    type(SUNMatrix), pointer :: B

    failure = 0

    ! protect A
    B => FSUNMatClone(A)

    failure = FSUNMatZero(B)
    if (failure /= 0) then
      call TEST_STATUS2(">>> FAILED test -- SUNMatZero returned  ", failure, myid)
      call FSUNMatDestroy(B)
      return
    end if

    ! A data should be a vector of zeros
    failure = check_matrix_entry(B, ZERO, tol);
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatZero check ", myid)
      call FSUNMatDestroy(B)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatZero ", myid)
    end if

    call FSUNMatDestroy(B)

  end function Test_FSUNMatZero

  integer(C_INT) function Test_FSUNMatCopy(A, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    real(C_DOUBLE)  :: tol = 10*UNIT_ROUNDOFF
    type(SUNMatrix) :: A
    type(SUNMatrix), pointer :: B

    failure = 0

    B => FSUNMatClone(A)

    ! copy matrix data
    failure = FSUNMatCopy(A, B)
    if (failure /= 0) then
      call TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned  ", failure, myid)
      call FSUNMatDestroy(B)
      failure = 1
      return
    end if

    ! check matrix entries
    failure = check_matrix(B, A, tol)
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatCopy ", myid)
      call FSUNMatDestroy(B)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatCopy ", myid)
    end if

    call FSUNMatDestroy(B)

  end function Test_FSUNMatCopy

  integer(C_INT) function Test_FSUNMatScaleAdd(A, I, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    real(C_DOUBLE)  :: tol = 10*UNIT_ROUNDOFF
    type(SUNMatrix) :: A, I
    type(SUNMatrix), pointer :: B

    failure = 0

    !
    ! Case 1: same sparsity/bandwith pattern
    !

    ! protect A
    B => FSUNMatClone(A)
    failure = FSUNMatCopy(A, B)
    if (failure /= 0) then
      call TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned  ", failure, myid)
      call FSUNMatDestroy(B)
      return
    end if

    ! fill vector data
    failure = FSUNMatScaleAdd(NEG_ONE, B, B)
    if (failure /= 0) then
      call TEST_STATUS2(">>> FAILED test -- SUNMatScaleAdd returned  ", failure, myid)
      call FSUNMatDestroy(B)
      return
    end if

    ! check matrix entries
    failure = check_matrix_entry(B, ZERO, tol)
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatScaleAdd ", myid)
      call FSUNMatDestroy(B)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatScaleAdd ", myid)
    end if

    call FSUNMatDestroy(B)

  end function Test_FSUNMatScaleAdd

  integer(C_INT) function Test_FSUNMatScaleAddI(A, I, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    real(C_DOUBLE)  :: tol = 10*UNIT_ROUNDOFF
    type(SUNMatrix) :: A, I
    type(SUNMatrix), pointer :: B

    failure = 0

    ! protect A
    B => FSUNMatClone(A)
    failure = FSUNMatCopy(I, B)
    if (failure /= 0) then
      call TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned  ", failure, myid)
      call FSUNMatDestroy(B)
      return
    end if

    ! fill vector data
    failure = FSUNMatScaleAddI(NEG_ONE, B)
    if (failure /= 0) then
      call TEST_STATUS2(">>> FAILED test -- SUNMatScaleAddI returned  ", failure, myid)
      call FSUNMatDestroy(B)
      return
    end if

    ! check matrix
    failure = check_matrix_entry(B, ZERO, tol)
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatScaleAddI check ", myid)
      call FSUNMatDestroy(B)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatScaleAddI ", myid)
    end if

    call FSUNMatDestroy(B)

  end function Test_FSUNMatScaleAddI

  integer(C_INT) function Test_FSUNMatMatvecSetup(A, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    type(SUNMatrix) :: A
    type(SUNMatrix_Ops), pointer :: ops

    failure = 0

    call c_f_pointer(A%ops, ops)

    if (c_associated(ops%matvecsetup)) then
      call TEST_STATUS("    PASSED test -- SUNMatMatvecSetup not implemented", myid)
      return
    end if

    failure = FSUNMatMatvecSetup(A)
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatMatvecSetup ", myid)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatMatvecSetup", myid)
    end if

  end function Test_FSUNMatMatvecSetup

  integer(C_INT) function Test_FSUNMatMatvec(A, x, y, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_nvector_mod
    use fsundials_matrix_mod

    implicit none

    type(SUNMatrix)              :: A
    type(SUNMatrix), pointer     :: B
    type(N_Vector)               :: x, y
    type(N_Vector),  pointer     :: z, w
    integer(C_INT)               :: myid
    real(C_DOUBLE)               :: tol = 100*UNIT_ROUNDOFF
    type(SUNMatrix_Ops), pointer :: ops

    failure = 0

    ! harder tests for square matrices
    if (is_square(A)) then

      ! protect A
      B => FSUNMatClone(A)
      failure = FSUNMatCopy(A, B)
      if (failure /= 0) then
        call TEST_STATUS2(">>> FAILED test -- SUNMatCopy returned  ", failure, myid)
        call FSUNMatDestroy(B)
        return
      end if

      ! compute matrix vector product
      failure = FSUNMatScaleAddI(THREE, B)
      if (failure /= 0) then
        call TEST_STATUS2(">>> FAILED test -- SUNMatScaleAddI returned  ", failure, myid)
        call FSUNMatDestroy(B)
        return
      end if

      z => FN_VClone(y)
      w => FN_VClone(y)

      ! Call the Setup function before the Matvec if it exists
      call c_f_pointer(A%ops, ops)
      if (c_associated(ops%matvecsetup)) then
        failure = FSUNMatMatvecSetup(B)
        if (failure /= 0) then
          call TEST_STATUS2(">>> FAILED test -- SUNMatMatvecSetup returned  ", failure, myid)
          call FSUNMatDestroy(B)
          return
        end if
      end if

      failure = FSUNMatMatvec(B,x,z)
      if (failure /= 0) then
        call TEST_STATUS2(">>> FAILED test -- SUNMatMatvec returned  ", failure, myid)
        call FSUNMatDestroy(B)
        return
      end if

      call FN_VLinearSum(THREE,y,ONE,x,w)

      failure = check_vector(w,z,tol)

      call FSUNMatDestroy(B)
      call FN_VDestroy(z)
      call FN_VDestroy(w)

    else

      z => FN_VClone(y)

      failure = FSUNMatMatvec(A,x,z)
      if (failure /= 0) then
        call TEST_STATUS2(">>> FAILED test -- SUNMatMatvec returned  ", failure, myid)
        return
      end if

      failure = check_vector(y,z,tol)
      call FN_VDestroy(z)

    end if

    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatMatvec check ", myid)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatMatvec ", myid)
    end if

  end function Test_FSUNMatMatvec

  integer(C_INT) function Test_FSUNMatSpace(A, myid) result(failure)
    use, intrinsic :: iso_c_binding
    use test_utilities
    use fsundials_matrix_mod

    implicit none

    integer(C_INT)  :: myid
    type(SUNMatrix) :: A
    integer(C_LONG) :: lenrw(1), leniw(1)

    failure = 0

    failure = FSUNMatSpace(A, lenrw, leniw);
    if (failure /= 0) then
      call TEST_STATUS(">>> FAILED test -- SUNMatSpace ", myid)
      return
    else
      call TEST_STATUS("    PASSED test -- SUNMatSpace ", myid)
    end if

  end function Test_FSUNMatSpace

end module test_sunmatrix
