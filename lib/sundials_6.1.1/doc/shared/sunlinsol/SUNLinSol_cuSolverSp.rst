..
   Programmer(s): Cody J. Balos @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.cuSolverSp:

The SUNLinSol_cuSolverSp_batchQR Module
=======================================

The SUNLinSol_cuSolverSp_batchQR implementation of the ``SUNLinearSolver`` class
is designed to be used with the SUNMATRIX_CUSPARSE matrix, and the NVECTOR_CUDA vector.
The header file to include when using this module is
``sunlinsol/sunlinsol_cusolversp_batchqr.h``. The installed library to link to
is ``libsundials_sunlinsolcusolversp.lib`` where ``.lib`` is typically
``.so`` for shared libraries and ``.a`` for static libraries.

.. warning::

   The SUNLinearSolver_cuSolverSp_batchQR module is experimental and subject to change.


.. _SUNLinSol.cuSolverSp.description:

SUNLinSol_cuSolverSp_batchQR description
----------------------------------------

The SUNLinearSolver_cuSolverSp_batchQR implementation provides an interface to
the batched sparse QR factorization method provided by the NVIDIA cuSOLVER library :cite:p:`cuSOLVER_site`.
The module is designed for solving block diagonal linear systems of the form

.. math::

   \begin{bmatrix}
      \mathbf{A_1} & 0 & \cdots & 0\\
      0 & \mathbf{A_2} & \cdots & 0\\
      \vdots & \vdots & \ddots & \vdots\\
      0 & 0 & \cdots & \mathbf{A_n}\\
   \end{bmatrix}
   x_j
   =
   b_j

where all block matrices :math:`\mathbf{A_j}` share the same sparsity pattern. The matrix
must be the ``SUNMatrix.cuSparse``.


.. _SUNLinSol.cuSolverSp.functions:

SUNLinSol_cuSolverSp_batchQR functions
--------------------------------------

The ``SUNLinearSolver_cuSolverSp_batchQR`` module defines implementations of
all "direct" linear solver operations listed in :numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_cuSolverSp_batchQR``

* ``SUNLinSolInitialize_cuSolverSp_batchQR`` -- this sets the
  ``first_factorize`` flag to 1

* ``SUNLinSolSetup_cuSolverSp_batchQR`` -- this always copies the
  relevant SUNMATRIX_SPARSE data to the GPU; if this is the first setup
  it will perform symbolic analysis on the system

* ``SUNLinSolSolve_cuSolverSp_batchQR`` -- this calls the
  ``cusolverSpXcsrqrsvBatched`` routine to perform factorization

* ``SUNLinSolLastFlag_cuSolverSp_batchQR``

* ``SUNLinSolFree_cuSolverSp_batchQR``


In addition, the module provides the following user-callable routines:

.. c:function:: SUNLinearSolver SUNLinSol_cuSolverSp_batchQR(N_Vector y, SUNMatrix A, cusolverHandle_t cusol, SUNContext sunctx)

   The function ``SUNLinSol_cuSolverSp_batchQR`` creates and allocates
   memory for a SUNLinearSolver object.

   **Arguments:**
      * *y* -- a vector for checking compatibility with the solver.
      * *A* -- a SUNMATRIX_cuSparse matrix for checking compatibility with the
        solver.
      * *cusol* -- cuSolverSp object to use.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNLinearSolver`` object. If either *A* or *y* are
      incompatible then this routine will return ``NULL``.

   **Notes:**
      This routine will perform consistency checks to ensure that it is
      called with consistent ``N_Vector``  and ``SUNMatrix``  implementations.
      These are currently limited to the SUNMATRIX_CUSPARSE matrix type
      and the NVECTOR_CUDA vector type. Since the SUNMATRIX_CUSPARSE matrix
      type is only compatible with the NVECTOR_CUDA the restriction is also
      in place for the linear solver. As additional compatible matrix and
      vector implementations are added to SUNDIALS, these will be included
      within this compatibility check.


.. c:function:: void SUNLinSol_cuSolverSp_batchQR_GetDescription(SUNLinearSolver LS, char **desc)

   The function ``SUNLinSol_cuSolverSp_batchQR_GetDescription``
   accesses the string description of the object (empty by default).


.. c:function:: void SUNLinSol_cuSolverSp_batchQR_SetDescription(SUNLinearSolver LS, const char *desc)

   The function ``SUNLinSol_cuSolverSp_batchQR_SetDescription``
   sets the string description of the object (empty by default).


.. c:function:: void SUNLinSol_cuSolverSp_batchQR_GetDeviceSpace(SUNLinearSolver S, size_t* cuSolverInternal, size_t* cuSolverWorkspace)

   The function ``SUNLinSol_cuSolverSp_batchQR_GetDeviceSpace``
   returns the cuSOLVER batch QR method internal buffer size, in bytes,
   in the argument ``cuSolverInternal`` and the cuSOLVER
   batch QR workspace buffer size, in bytes, in the agrument
   ``cuSolverWorkspace``. The size of the internal buffer is
   proportional to the number of matrix blocks while the size
   of the workspace is almost independent of the number of blocks.


.. _SUNLinSol.cuSolverSp.content:

SUNLinSol_cuSolverSp_batchQR content
------------------------------------

The SUNLinSol_cuSolverSp_batchQR module defines the *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_cuSolverSp_batchQR {
      int                last_flag;       /* last return flag                          */
      booleantype        first_factorize; /* is this the first factorization?          */
      size_t             internal_size;   /* size of cusolver buffer for Q and R       */
      size_t             workspace_size;  /* size of cusolver memory for factorization */
      cusolverSpHandle_t cusolver_handle; /* cuSolverSp context                        */
      csrqrInfo_t        info;            /* opaque cusolver data structure            */
      void*              workspace;       /* memory block used by cusolver             */
      const char*        desc;            /* description of this linear solver         */
   };
