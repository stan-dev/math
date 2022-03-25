.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _NVectors.ARKODE:

NVECTOR functions required by ARKODE
==========================================

In :numref:`NVectors.ARKODE.Table` below, we list the vector functions in
the ``N_Vector`` module that are called within the ARKODE package.  The
table also shows, for each function, which ARKODE module uses the function.
The ARKSTEP and ERKSTEP columns show function usage within the main
time-stepping modules and the shared ARKODE infrastructure,  while the
remaining columns show function usage within the ARKLS linear solver
interface, the ARKBANDPRE and ARKBBDPRE preconditioner modules.

Note that for ARKLS we only list the ``N_Vector`` routines used
directly by ARKLS, each ``SUNLinearSolver`` module may have additional
requirements that are not listed here.  In addition, specific
``SUNNonlinearSolver`` modules attached to ARKODE may have additional
``N_Vector`` requirements.  For additional requirements by specific
``SUNLinearSolver`` and ``SUNNonlinearSolver`` modules, please see the
accompanying sections :numref:`SUNLinSol` and :numref:`SUNNonlinSol`.

At this point, we should emphasize that the user does not need to know
anything about ARKODE's usage of vector functions in order to use
ARKODE.  Instead, this information is provided primarily for users
interested in constructing a custom ``N_Vector`` module.  We note that
a number of ``N_Vector`` functions from the section
:numref:`NVectors.Description` are not listed in the above table.
Therefore a user-supplied ``N_Vector`` module for ARKODE could safely
omit these functions from their implementation (although
some may be needed by ``SUNNonlinearSolver`` or ``SUNLinearSolver``
modules).


.. _NVectors.ARKODE.Table:
.. table:: List of vector functions usage by ARKODE code modules

   ========================================  =======  =======  =======  =====  ==========  =========
   Routine                                   ARKSTEP  ERKSTEP  MRISTEP  ARKLS  ARKBANDPRE  ARKBBDPRE
   ========================================  =======  =======  =======  =====  ==========  =========
   :c:func:`N_VGetLength`                                               4
   :c:func:`N_VAbs`                          X        X
   :c:func:`N_VAddConst`                     X        X
   :c:func:`N_VClone`                        X        X        X        X
   :c:func:`N_VCloneEmpty`
   :c:func:`N_VConst`                        X        X        X        X
   :c:func:`N_VDestroy`                      X        X        X        X
   :c:func:`N_VDiv`                          X        X
   :c:func:`N_VGetArrayPointer`                                         1      X           X
   :c:func:`N_VInv`                          X        X
   :c:func:`N_VLinearSum`                    X        X        X        X
   :c:func:`N_VMaxNorm`                      X        X
   :c:func:`N_VMin`                          X        X
   :c:func:`N_VScale`                        X        X        X        X      X           X
   :c:func:`N_VSetArrayPointer`                                         1
   :c:func:`N_VSpace`\ :sup:`2`              X        X        X        X      X           X
   :c:func:`N_VWrmsNorm`                     X        X        X        X      X           X
   :c:func:`N_VLinearCombination`\ :sup:`3`  X        X        X
   :c:func:`N_VMinQuotient`\ :sup:`5`        X        X
   :c:func:`N_VConstrMask`\ :sup:`5`         X        X
   :c:func:`N_VCompare`\ :sup:`5`            X        X
   ========================================  =======  =======  =======  =====  ==========  =========

Special cases (numbers match markings in table):

1. This is only required with the :ref:`SUNMATRIX_DENSE <SUNMatrix.Dense>` or
   :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` modules,
   where the default difference-quotient Jacobian approximation is used.

2. The :c:func:`N_VSpace()` function is only informational, and will
   only be called if provided by the ``N_Vector`` implementation.

3. The :c:func:`N_VLinearCombination()` function is in fact optional;
   if it is not supplied then :c:func:`N_VLinearSum()` will be used instead.

4. The :c:func:`N_VGetLength()` function is only required when an iterative or
   matrix iterative ``SUNLinearSolver`` module is used.

5. The functions :c:func:`N_VMinQuotient`, :c:func:`N_VConstrMask`, and
   :c:func:`N_VCompare` are only used when inequality constraints are enabled
   and may be omitted if this feature is not used.
