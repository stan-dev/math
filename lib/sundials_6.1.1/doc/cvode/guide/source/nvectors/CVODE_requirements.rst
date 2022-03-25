.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.CVODE:

NVECTOR functions used by CVODE
===============================

In :numref:`NVectors.CVODE.Table` below, we list the vector functions in the
NVECTOR module used within the CVODE package.
The table also shows, for each function, which of the code modules uses
the function. The CVODE column shows function usage within the main
integrator module, while the remaining columns show function usage
within each of the CVODE linear solver interfaces, the CVBANDPRE and
CVBBDPRE preconditioner modules. Here
CVLS stands for the generic linear solver interface in CVODE,
and CVDIAG stands for the diagonal linear solver interface in CVODE.

At this point, we should emphasize that the CVODE user does not need to know
anything about the usage of vector functions by the CVODE code modules in order
to use CVODE. The information is presented as an implementation detail for the
interested reader.

.. _NVectors.CVODE.Table:
.. table:: List of vector functions usage by CVODE code modules

   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   |                                | CVODE       | CVLS       | CVDIAG       | CVBANDPRE       | CVBBDPRE       |
   +================================+=============+============+==============+=================+================+
   | :c:func:`N_VGetVectorID`       |             |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VGetLength`         |             | 4          |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VClone`             | x           | x          | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VCloneEmpty`        |             | 1          |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VDestroy`           | x           | x          | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VSpace`             | x           | 2          |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VGetArrayPointer`   |             | 1          |              | x               | x              |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VSetArrayPointer`   |             | 1          |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VLinearSum`         | x           | x          | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VConst`             | x           | x          |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VProd`              | x           |            | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VDiv`               | x           |            | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VScale`             | x           | x          | x            | x               | x              |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VAbs`               | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VInv`               | x           |            | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VAddConst`          | x           |            | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VMaxNorm`           | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VWrmsNorm`          | x           | x          |              | x               | x              |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VMin`               | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VMinQuotient`       | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VConstrMask`        | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VCompare`           | x           |            | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VInvTest`           |             |            | x            |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VLinearCombination` | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VScaleAddMulti`     | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VDotProdMulti`      | 3           | 3          |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+
   | :c:func:`N_VScaleVectorArray`  | x           |            |              |                 |                |
   +--------------------------------+-------------+------------+--------------+-----------------+----------------+


Special cases (numbers match markings in table):

1. These routines are only required if an internal
   difference-quotient routine for constructing :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>`
   or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices is used.

2. This routine is optional, and is only used in estimating
   space requirements for CVODE modules for user feedback.

3. The optional function ``N_VDotProdMulti`` is only used in the
   SUNNONLINSOL_FIXEDPOINT module, or when Classical Gram-Schmidt is
   enabled with SPGMR or SPFGMR. The remaining operations from
   :numref:`NVectors.Ops.Fused` -- :numref:`NVectors.Ops.Array` not listed above
   are unused and a user-supplied NVECTOR module for CVODE could
   omit these operations.

4. This routine is only used when an iterative or matrix iterative
   SUNLINSOL module is supplied to CVODE.

Each SUNLINSOL object may require additional NVECTOR routines
not listed in the table above. Please see the the relevant
descriptions of these modules in :numref:`SUNLinSol` for
additional detail on their NVECTOR requirements.

The vector functions listed in :numref:`NVectors.Ops` that are *not* used by
CVODE are: :c:func:`N_VWL2Norm`, :c:func:`N_VDotProd`, :c:func:`N_VL1Norm`,
:c:func:`N_VWrmsNormMask`, and :c:func:`N_VGetCommunicator`. Therefore, a
user-supplied NVECTOR module for CVODE could omit these functions (although
some may be needed by SUNNONLINSOL or SUNLINSOL modules). The functions
:c:func:`N_MinQuotient`, :c:func:`N_VConstrMask`, and :c:func:`N_VCompare`
are only used when constraint checking is enabled and may be omitted if this
feature is not used.
