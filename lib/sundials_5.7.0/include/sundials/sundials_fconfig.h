!
! ----------------------------------------------------------------- 
! Programmer(s): Daniel R. Reynolds @ SMU
!-----------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2021, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ------------------------------------------------------------------
! SUNDIALS fortran configuration input
! ------------------------------------------------------------------

!     Define precision of SUNDIALS data type 'realtype' as Fortran 
!     parameter "REALTYPE"
!
!     Depending on the precision level, this value will be one of
!          4  (SUNDIALS_SINGLE_PRECISION)
!          8  (SUNDIALS_DOUBLE_PRECISION)
!         16  (SUNDIALS_EXTENDED_PRECISION)
!
integer REALTYPE
parameter (REALTYPE=8)

!     Define type of vector indices in SUNDIALS 'sunindextype' as 
!     the Fortran parameter "SUNINDEXTYPE"
!
!     Depending on the user choice of indextype, this will be one of
!          4  (32BIT)
!          8  (64BIT)
!
integer SUNINDEXTYPE
parameter (SUNINDEXTYPE=8)

!     If building with MPI enabled, define the logical flag 
!     "SUNDIALS_MPI_COMM_F2C" indicating whether the user can specify
!     a different MPI communicator than MPI_COMM_WORLD to FNVInitP
!
!          .true.   (communicator can differ from MPI_COMM_WORLD)
!          .false.  (communicator must be MPI_COMM_WORLD)
!
logical SUNDIALS_MPI_COMM_F2C
parameter (SUNDIALS_MPI_COMM_F2C=.true.)
