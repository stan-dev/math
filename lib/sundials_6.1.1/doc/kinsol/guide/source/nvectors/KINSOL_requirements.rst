.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.KINSOL:

NVECTOR functions used by KINSOL
================================

In :numref:`NVectors.KINSOL.Table` below, we list the vector functions used in the ``N_Vector`` module used
by the KINSOL package. The table also shows, for each function, which of the code modules uses the
function. The KINSOL column shows function usage within the main integrator module, while the remaining
columns show function usage within the KINLS linear solvers interface, and the KINBBDPRE preconditioner
module.

At this point, we should emphasize that the KINSOL user does not need to know anything about the usage
of vector functions by the KINSOL code modules in order to use KINSOL. The information is presented as an
implementation detail for the interested reader.

.. _NVectors.KINSOL.Table:
.. table:: List of vector functions usage by KINSOL code modules

   ==============================  ======  =====  =========
   Function name                   KINSOL  KINLS  KINBBDPRE
   ==============================  ======  =====  =========
   :c:func:`N_VGetVectorID`
   :c:func:`N_VGetLength`          4
   :c:func:`N_VClone`              x       x
   :c:func:`N_VCloneEmpty`
   :c:func:`N_VDestroy`            x       x
   :c:func:`N_VSpace`              x       2
   :c:func:`N_VGetArrayPointer`    1       x
   :c:func:`N_VSetArrayPointer`    1
   :c:func:`N_VLinearSum`          x       x
   :c:func:`N_VConst`              x
   :c:func:`N_VProd`               x       x
   :c:func:`N_VDiv`                x
   :c:func:`N_VScale`              x       x      x
   :c:func:`N_VAbs`                x
   :c:func:`N_VInv`                x
   :c:func:`N_VDotProd`            x       x
   :c:func:`N_VMaxNorm`            x
   :c:func:`N_VMin`                x
   :c:func:`N_VWL2Norm`            x       x
   :c:func:`N_VL1Norm`             3
   :c:func:`N_VConstrMask`         x
   :c:func:`N_VMinQuotient`        x
   :c:func:`N_VLinearCombination`  x       x
   :c:func:`N_VDotProdMulti`       x
   ==============================  ======  =====  =========

Special cases (numbers match markings in table):

#. These routines are only required if an internal difference-quotient routine
   for constructing :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
   :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices is used.

#. This routine is optional, and is only used in estimating space requirements
   for IDA modules for user feedback.

#. These routines are only required if the internal difference-quotient routine
   for approximating the Jacobian-vector product is used.

#. This routine is only used when an iterative ``SUNLinearSolver`` module that
   does not support the :c:func:`SUNLinSolSetScalingVectors` routine is supplied to KINSOL.

Each ``SUNLinearSolver`` object may require additional ``N_Vector`` routines not
listed in the table above. Please see the the relevant descriptions of these
modules in :numref:`SUNLinSol` for additional detail on their ``N_Vector``
requirements.

The vector functions listed in :numref:`NVectors.Ops` that are *not* used by
KINSOL are :c:func:`N_VAddConst`, :c:func:`N_VWrmsNorm`, :c:func:`N_VWrmsNormMask`,
:c:func:`N_VCompare`, :c:func:`N_VInvTest`, and :c:func:`N_VGetCommunicator`. Therefore a
user-supplied ``N_Vector`` module for KINSOL could omit these functions.

The optional function :c:func:`N_VLinearCombination` is only used when Anderson
acceleration is enabled or the SPBCG, SPTFQMR, SPGMR, or SPFGMR linear solvers
are used. :c:func:`N_VDotProd` is only used when Anderson acceleration is enabled or
Classical Gram-Schmidt is used with SPGMR or SPFGMR. The remaining operations
from :numref:`NVectors.Ops.Fused` and :numref:`NVectors.Ops.Array` are unused and a
user-supplied ``N_Vector`` module for KINSOL could omit these operations.
