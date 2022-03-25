..
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNLinSol.MagmaDense:

The SUNLinSol_MagmaDense Module
======================================

The SUNLinearSolver_MagmaDense implementation of the ``SUNLinearSolver`` class is
designed to be used with the SUNMATRIX_MAGMADENSE matrix, and a GPU-enabled
vector. The header file to include when using this module is
``sunlinsol/sunlinsol_magmadense.h``. The installed library to link to is
``libsundials_sunlinsolmagmadense.lib`` where ``lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

.. warning::

   The SUNLinearSolver_MagmaDense module is experimental and subject to change.


SUNLinearSolver_MagmaDense Description
---------------------------------------

The SUNLinearSolver_MagmaDense implementation provides an interface to the dense
LU and dense batched LU methods in the `MAGMA <https://icl.utk.edu/magma/index.html>`_
linear algebra library :cite:p:`magma_ref`. The batched LU methods are leveraged when solving block
diagonal linear systems of the form

.. math::

   \begin{bmatrix}
     \mathbf{A_0} & 0 & \cdots & 0\\
     0 & \mathbf{A_1} & \cdots & 0\\
     \vdots & \vdots & \ddots & \vdots\\
     0 & 0 & \cdots & \mathbf{A_{n-1}}\\
   \end{bmatrix}
   x_j
   =
   b_j.


SUNLinearSolver_MagmaDense Functions
-------------------------------------

The SUNLinearSolver_MagmaDense module defines implementations of all "direct"
linear solver operations listed in :numref:`SUNLinSol.API`:

* ``SUNLinSolGetType_MagmaDense``
* ``SUNLinSolInitialize_MagmaDense``
* ``SUNLinSolSetup_MagmaDense``
* ``SUNLinSolSolve_MagmaDense``
* ``SUNLinSolLastFlag_MagmaDense``
* ``SUNLinSolFree_MagmaDense``

In addition, the module provides the following user-callable routines:


.. c:function:: SUNLinearSolver SUNLinSol_MagmaDense(N_Vector y, SUNMatrix A, SUNContext sunctx)

   This constructor function creates and allocates memory for a
   ``SUNLinearSolver`` object.

   **Arguments:**
      * *y* -- a vector for checking compatibility with the solver.
      * *A* -- a SUNMATRIX_MAGMADENSE matrix for checking compatibility with the
        solver.
      * *sunctx* -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)

   **Return value:**
      If successful, a ``SUNLinearSolver`` object. If either *A* or *y* are
      incompatible then this routine will return ``NULL``. This routine analyzes
      the input matrix and vector to determine the linear system size and to
      assess compatibility with the solver.

.. c:function:: int SUNLinSol_MagmaDense_SetAsync(SUNLinearSolver LS, booleantype onoff)

   This function can be used to toggle the linear solver between asynchronous
   and synchronous modes. In asynchronous mode (default), SUNLinearSolver
   operations are asynchronous with respect to the host. In synchronous mode,
   the host and GPU device are synchronized prior to the operation returning.

   **Arguments:**
      * *LS* -- a SUNLinSol_MagmaDense object
      * *onoff* -- 0 for synchronous mode or 1 for asynchronous mode (default 1)

   **Return value:**
      * ``SUNLS_SUCCESS`` if successful
      * ``SUNLS_MEM_NULL`` if *LS* is ``NULL``


SUNLinearSolver_MagmaDense Content
-----------------------------------

The SUNLinearSolver_MagmaDense module defines the object *content* field of a
``SUNLinearSolver`` to be the following structure:

.. code-block:: c

   struct _SUNLinearSolverContent_MagmaDense {
     int             last_flag;
     booleantype     async;
     sunindextype    N;
     SUNMemory       pivots;
     SUNMemory       pivotsarr;
     SUNMemory       dpivotsarr;
     SUNMemory       infoarr;
     SUNMemory       rhsarr;
     SUNMemoryHelper memhelp;
     magma_queue_t   q;
   };
