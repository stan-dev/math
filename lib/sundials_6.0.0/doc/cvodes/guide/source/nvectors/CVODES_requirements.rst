.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.CVODES:

NVECTOR functions used by CVODES
================================

In :numref:`NVectors.CVODES.Table` below, we list the vector functions in the
``N_Vector`` module used within the CVODES package. The table also shows, for
each function, which of the code modules uses the function. The CVODES column
shows function usage within the main integrator module, while the remaining
columns show function usage within each of the CVODES linear solver interfaces,
the CVBANDPRE and CVBBDPRE preconditioner modules, and the CVODES adjoint
sensitivity module (denoted here by CVODEA). Here CVLS stands for the
generic linear solver interface in CVODES, and CVDIAG stands for the diagonal
linear solver interface in CVODES.

At this point, we should emphasize that the CVODES user does not need to know
anything about the usage of vector functions by the CVODES code modules in order
to use CVODES. The information is presented as an implementation detail for the
interested reader.

.. _NVectors.CVODES.Table:
.. table:: List of vector functions usage by CVODES code modules

  ========================================= ====== ==== ====== ========= ======== ======
  \                                         CVODES CVLS CVDIAG CVBANDPRE CVBBDPRE CVODEA
  ========================================= ====== ==== ====== ========= ======== ======
  :c:func:`N_VGetVectorID`
  :c:func:`N_VGetLength`                             4
  :c:func:`N_VClone`                          x      x    x                         x
  :c:func:`N_VCloneEmpty`                            1
  :c:func:`N_VDestroy`                        x      x    x                         x
  :c:func:`N_VCloneVectorArray`               x                                     x
  :c:func:`N_VDestroyVectorArray`             x                                     x
  :c:func:`N_VSpace`                          x      2
  :c:func:`N_VGetArrayPointer`                       1           x         x
  :c:func:`N_VSetArrayPointer`                       1
  :c:func:`N_VLinearSum`                      x      x    x                         x
  :c:func:`N_VConst`                          x      x
  :c:func:`N_VProd`                           x           x
  :c:func:`N_VDiv`                            x           x
  :c:func:`N_VScale`                          x      x    x      x         x        x
  :c:func:`N_VAbs`                            x
  :c:func:`N_VInv`                            x           x
  :c:func:`N_VAddConst`                       x           x
  :c:func:`N_VMaxNorm`                        x
  :c:func:`N_VWrmsNorm`                       x      x           x         x
  :c:func:`N_VMin`                            x
  :c:func:`N_MinQuotient`                     x
  :c:func:`N_VConstrMask`                     x
  :c:func:`N_VCompare`                        x           x
  :c:func:`N_VInvTest`                                    x
  :c:func:`N_VLinearCombination`              x
  :c:func:`N_VScaleAddMulti`                  x
  :c:func:`N_VDotProdMulti`                   3      3
  :c:func:`N_VLinearSumVectorArray`           x
  :c:func:`N_VScaleVectorArray`               x
  :c:func:`N_VConstVectorArray`               x
  :c:func:`N_VWrmsNormVectorArray`            x
  :c:func:`N_VScaleAddMultiVectorArray`       x
  :c:func:`N_VLinearCombinationVectorArray`   x
  ========================================= ====== ==== ====== ========= ======== ======


Special cases (numbers match markings in table):

1. These routines are only required if an internal difference-quotient routine
   for constructing :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
   :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` Jacobian matrices is used.

2. This routine is optional, and is only used in estimating space requirements
   for CVODES modules for user feedback.

3. The optional function :c:func:`N_VDotProdMulti` is only used in the
   ``SUNNONLINSOL_FIXEDPOINT`` module, or when Classical
   Gram-Schmidt is enabled with SPGMR or SPFGMR.

4. This routine is only used when an iterative or matrix iterative
   ``SUNLinearSolver`` module is supplied to CVODES.

Each ``SUNLinearSolver`` object may require additional ``N_Vector`` routines not
listed in the table above. Please see the the relevant descriptions of these
modules in :numref:`SUNLinSol` for additional detail on their ``N_Vector``
requirements.

The remaining operations from :numref:`NVectors.Ops` not listed above are unused
and a user-supplied ``N_Vector`` module for CVODES could omit these operations
(although some may be needed by ``SUNNonlinearSolver`` or ``SUNLinearSolver``
modules). The functions :c:func:`N_VMinQuotient`, :c:func:`N_VConstrMask`, and
:c:func:`N_VCompare` are only used when constraint checking is enabled and may
be omitted if this feature is not used.
