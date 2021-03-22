! -----------------------------------------------------------------
! Programmer(s): Cody J. Balos @ LLNL
! -----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
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

contains

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
