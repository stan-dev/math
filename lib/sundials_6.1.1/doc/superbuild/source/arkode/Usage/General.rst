.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
                  Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Usage.Headers:

Access to library and header files
===========================================

At this point, it is assumed that the installation of ARKODE,
following the procedure described in :numref:`Installation`,
has been completed successfully.

Regardless of where the user's application program resides, its
associated compilation and load commands must make reference to the
appropriate locations for the library and header files required by
ARKODE. The relevant library files are

- ``libdir/libsundials_arkode.lib``,
- ``libdir/libsundials_nvec*.lib``,

where the file extension ``.lib`` is typically ``.so`` for shared
libraries and ``.a`` for static libraries.  The relevant header files
are located in the subdirectories

- ``incdir/include/arkode``
- ``incdir/include/sundials``
- ``incdir/include/nvector``
- ``incdir/include/sunmatrix``
- ``incdir/include/sunlinsol``
- ``incdir/include/sunnonlinsol``

The directories ``libdir`` and ``incdir`` are the installation library
and include directories, respectively.  For a default installation,
these are ``instdir/lib`` and ``instdir/include``, respectively, where
``instdir`` is the directory where SUNDIALS was installed (see
:numref:`Installation` for further details).


When using ARKODE, the calling program must include several header
files so that various macros and data types can be used. One of the
following header files is always required:

- ``arkode/arkode_arkstep.h``, the main header file for the ARKStep
  time-stepping module.

- ``arkode/arkode_erkstep.h``, the main header file for the ERKStep
  time-stepping module.

- ``arkode/arkode_mristep.h``, the main header file for the MRIStep
  time-stepping module.

Each of these define several types and various constants, include
function prototypes, and include the shared ``arkode/arkode.h`` and
``arkode/arkode_ls.h`` header files.

Note that ``arkode.h`` includes ``sundials_types.h`` directly, which
defines the types ``realtype``,  ``sunindextype``, and ``booleantype``
and the constants ``SUNFALSE`` and ``SUNTRUE``, so a user program does
not need to include ``sundials_types.h`` directly.

Additionally, the calling program must also include an NVECTOR
implementation header file, of the form ``nvector/nvector_***.h``,
corresponding to the user's preferred data layout and form of
parallelism.  See :numref:`NVectors` for details for the
appropriate name.  This file in turn includes the header file
``sundials_nvector.h`` which defines the abstract ``N_Vector`` data
type.

If the user wishes to manually select between any of the pre-defined
ERK or DIRK Butcher tables (for ARKStep, ERKStep, or as the basis for
an MIS method), these are defined through a set of constants
that are enumerated in the header files ``arkode/arkode_butcher_erk.h``
and ``arkode/arkode_butcher_dirk.h``, or if a user wishes to manually
specify one or more Butcher tables, the corresponding ``ARKodeButcherTable``
structure is defined in ``arkode/arkode_butcher.h``. Alternatively,
for MRIStep, slow-to-fast coupling coefficient tables are enumerated in the
header file ``arkode/arkode_mristp.h``, or if a user wishes to manually specify
a coupling table, the corresponding ``MRIStepCouplingMem`` structure is defined
in ``arkode/arkode_mristep.h``.

If the user includes a non-trivial implicit component to their ODE
system in ARKStep, or if the slow time scale for MRIStep should be treated
implicitly, then each implicit stage will require a nonlinear solver for
the resulting system of algebraic equations -- the default for this is a
modified or inexact Newton iteration, depending on the user's choice of
linear solver.  If using a non-default nonlinear solver
module, or when interacting with a SUNNONLINSOL module directly, the
calling program must also include a SUNNONLINSOL header file, of the
form ``sunnonlinsol/sunnonlinsol_***.h`` where ``***`` is the name of
the nonlinear solver module (see :numref:`SUNNonlinSol` for
more information).  This file in turn includes the header file
``sundials_nonlinearsolver.h`` which defines the abstract
``SUNNonlinearSolver`` data type.

If using a nonlinear solver that requires the solution of a linear
system of the form :math:`\mathcal{A}x=b` (e.g., the default Newton
iteration), then a linear solver module header file will also be
required.  Similarly, if the ODE system in ARKStep involves a non-identity
mass matrix :math:`M \ne I`, then each time step will require a linear
solver for systems of the form :math:`Mx=b`.  The header files
corresponding to the SUNDIALS-provided linear solver modules available
for use with ARKODE are:

- Direct linear solvers:

  - ``sunlinsol/sunlinsol_dense.h``,
    which is used with the dense linear solver module,
    SUNLINSOL_DENSE;

  - ``sunlinsol/sunlinsol_band.h``,
    which is used with the banded linear solver module,
    SUNLINSOL_BAND;

  - ``sunlinsol/sunlinsol_lapackdense.h``,
    which is used with the LAPACK dense linear solver module,
    SUNLINSOL_LAPACKDENSE;

  - ``sunlinsol/sunlinsol_lapackband.h``,
    which is used with the LAPACK banded linear solver module,
    SUNLINSOL_LAPACKBAND;

  - ``sunlinsol/sunlinsol_klu.h``,
    which is used with the KLU sparse linear solver module,
    SUNLINSOL_KLU;

  - ``sunlinsol/sunlinsol_superlumt.h``,
    which is used with the SuperLU_MT sparse linear solver module,
    SUNLINSOL_SUPERLUMT;

  - ``sunlinsol/sunlinsol_superludist.h``,
    which is used with the SuperLU_DIST parallel sparse linear solver module,
    SUNLINSOL_SUPERLUDIST;

  - ``sunlinsol/sunlinsol_cusolversp_batchqr.h``,
    which is used with the batched sparse QR factorization method provided
    by the NVDIA cuSOLVER library, SUNLINSOL_CUSOLVERSP_BATCHQR;

- Iterative linear solvers:

  - ``sunlinsol/sunlinsol_spgmr.h``,
    which is used with the scaled, preconditioned GMRES Krylov linear
    solver module, SUNLINSOL_SPGMR;

  - ``sunlinsol/sunlinsol_spfgmr.h``,
    which is used with the scaled, preconditioned FGMRES Krylov linear
    solver module, SUNLINSOL_SPFGMR;

  - ``sunlinsol/sunlinsol_spbcgs.h``,
    which is used with the scaled, preconditioned Bi-CGStab Krylov
    linear solver module, SUNLINSOL_SPBCGS;

  - ``sunlinsol/sunlinsol_sptfqmr.h``,
    which is used with the scaled, preconditioned TFQMR Krylov linear
    solver module, SUNLINSOL_SPTFQMR;

  - ``sunlinsol/sunlinsol_pcg.h``,
    which is used with the scaled, preconditioned CG Krylov linear
    solver module, SUNLINSOL_PCG;

The header files for the SUNLINSOL_DENSE and SUNLINSOL_LAPACKDENSE
linear solver modules include the file
``sunmatrix/sunmatrix_dense.h``, which defines the SUNMATRIX_DENSE
matrix module, as well as various functions and macros for acting on
such matrices.

The header files for the SUNLINSOL_BAND and SUNLINSOL_LAPACKBAND
linear solver modules include the file ``sunmatrix/sunmatrix_band.h``,
which defines the SUNMATRIX_BAND matrix module, as well as various
functions and macros for acting on such matrices.

The header files for the SUNLINSOL_KLU and SUNLINSOL_SUPERLUMT linear
solver modules include the file ``sunmatrix/sunmatrix_sparse.h``,
which defines the SUNMATRIX_SPARSE matrix module, as well as various
functions and macros for acting on such matrices.

The header file for the SUNLINSOL_CUSOLVERSP_BATCHQR
linear solver module includes the file ``sunmatrix/sunmatrix_cusparse.h``,
which defines the SUNMATRIX_CUSPARSE matrix module, as well as various
functions for acting on such matrices.

The header file for the SUNLINSOL_SUPERLUDIST
linear solver module includes the file ``sunmatrix/sunmatrix_slunrloc.h``,
which defines the SUNMATRIX_SLUNRLOC matrix module, as well as various
functions for acting on such matrices.

The header files for the Krylov iterative solvers include the file
``sundials/sundials_iterative.h``, which enumerates the
preconditioning type and (for the SPGMR and SPFGMR solvers) the
choices for the Gram-Schmidt orthogonalization process.

Other headers may be needed, according to the choice of
preconditioner, etc.  For example, if preconditioning for an iterative
linear solver were performed using the ARKBBDPRE module, the header
``arkode/arkode_bbdpre.h`` is needed to access the preconditioner
initialization routines.


.. include:: ../../../../shared/Types.rst
