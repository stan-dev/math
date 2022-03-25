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
! A helper module for the fortran tests
! -----------------------------------------------------------------

module test_utilities

    use, intrinsic :: iso_c_binding
    use fsundials_context_mod
    implicit none

    real(C_DOUBLE), parameter :: UNIT_ROUNDOFF = epsilon(1.0d0)

    real(C_DOUBLE) :: NEG_TWO  = -2.0d0
    real(C_DOUBLE) :: NEG_ONE  = -1.0d0
    real(C_DOUBLE) :: NEG_HALF = -0.50d0
    real(C_DOUBLE) :: ZERO     = 0.0d0
    real(C_DOUBLE) :: HALF     = 0.5d0
    real(C_DOUBLE) :: ONE      = 1.0d0
    real(C_DOUBLE) :: TWO      = 2.0d0
    real(C_DOUBLE) :: THREE    = 3.0d0
    real(C_DOUBLE) :: FOUR     = 4.0d0
    real(C_DOUBLE) :: FIVE     = 5.0d0
    real(C_DOUBLE) :: SIX      = 6.0d0

    type(C_PTR)    :: sunctx

contains

  subroutine Test_Init(comm)
    implicit none
    type(C_PTR), value :: comm
    integer(C_INT)     :: retval

    retval = FSUNContext_Create(comm, sunctx)
    if (retval /= 0) then
      print *, 'ERROR in Test_Init: FSUNContext_Create returned nonzero'
      stop 1
    end if

  end subroutine

  subroutine Test_Finalize()
    implicit none
    integer(C_INT) :: retval

    retval = FSUNContext_Free(sunctx)

  end subroutine

  integer(C_INT) function FNEQTOL(a, b, tol) result(nequal)
    implicit none
    real(C_DOUBLE) :: a, b, tol

    if (a /= a) then
      nequal = 1
    else if ((abs(a-b)/abs(b)) > tol) then
      nequal = 1
    else
      nequal = 0
    end if

  end function FNEQTOL

  integer(C_INT) function FNEQ(a, b) result(nequal)
    implicit none
    real(C_DOUBLE) :: a, b

    if (a /= a) then
      nequal = 1
    else if ((abs(a-b)/abs(b)) > (10*UNIT_ROUNDOFF)) then
      nequal = 1
    else
      nequal = 0
    end if
  end function FNEQ

end module
