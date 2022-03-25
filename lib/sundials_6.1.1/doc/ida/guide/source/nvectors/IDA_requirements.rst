.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.IDA:

NVECTOR functions used by IDA
=============================

In :numref:`NVectors.IDA.Table` below, we list the vector functions used in the ``N_Vector`` module used
by the IDA package. The table also shows, for each function, which of the code modules uses the
function. The IDA column shows function usage within the main integrator module, while the remaining
columns show function usage within the IDALS linear solvers interface, and the IDABBDPRE preconditioner
module.

At this point, we should emphasize that the IDA user does not need to know anything about the usage
of vector functions by the IDA code modules in order to use IDA. The information is presented as an
implementation detail for the interested reader.

.. _NVectors.IDA.Table:
.. table:: List of vector functions usage by IDA code modules

   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   |                                   |    IDA                  |   IDALS                 |    IDABBDPRE            |
   +===================================+=========================+=========================+=========================+
   | :c:func:`N_VGetVectorID`          |                         |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VGetLength`            |                         | 4                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VClone`                | x                       | x                       | x                       |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VCloneEmpty`           |                         | 1                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VDestroy`              | x                       | x                       | x                       |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VSpace`                | x                       | 2                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VGetArrayPointer`      |                         | 1                       | x                       |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VSetArrayPointer`      |                         | 1                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VLinearSum`            | x                       | x                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VConst`                | x                       | x                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VProd`                 | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VDiv`                  | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VScale`                | x                       | x                       | x                       |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VAbs`                  | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VInv`                  | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VAddConst`             | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VMaxNorm`              | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VWrmsNorm`             | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VMin`                  | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VMinQuotient`          | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VConstrMask`           | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VWrmsNormMask`         | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VCompare`              | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VLinearCombination`    | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VScaleAddMulti`        | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VDotProdMulti`         |                         | 3                       |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VLinearSumVectorArray` | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+
   | :c:func:`N_VScaleVectorArray`     | x                       |                         |                         |
   +-----------------------------------+-------------------------+-------------------------+-------------------------+

Special cases (numbers match markings in table):

#. These routines are only required if an internal difference-quotient routine for constructing
   :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices is used.

#. This routine is optional, and is only used in estimating space requirements for IDA modules for
   user feedback.

#. The optional function ``N_VDotProdMulti`` is only used when Classical Gram-Schmidt is enabled
   with SPGMR or SPFGMR. The remaining operations from Tables :numref:`NVectors.Ops.Fused` and
   :numref:`NVectors.Ops.Array` not listed above are unused and a user-supplied
   ``N_Vector`` module for IDA could omit these operations.

#. This routine is only used when an iterative or matrix iterative ``SUNLinearSolver`` module is
   supplied to IDA.

Of the functions listed in :numref:`NVectors.Ops`, :c:func:`N_VWL2Norm`, :c:func:`N_VL1Norm`,
:c:func:`N_VInvTest`, and :c:func:`N_VGetCommunicator` are *not* used by IDA. Therefore a user-supplied
``N_Vector`` module for IDA could omit these functions (although some may be needed by
``SUNNonlinearSolver`` or ``SUNLinearSolver`` modules).
