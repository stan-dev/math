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

.. _NVectors.MPIPlusX:

The NVECTOR_MPIPLUSX Module
===========================

The NVECTOR_MPIPLUSX module is designed to facilitate the MPI+X
paradigm, where X is some form of on-node (local) parallelism
(e.g. OpenMP, CUDA).  This paradigm is becoming increasingly popular
with the rise of heterogeneous computing architectures.

The NVECTOR_MPIPLUSX implementation is designed to work with any
NVECTOR that implements the minimum "standard" set
of operations in :numref:`NVectors.Ops.Standard`.
However, it is not recommended to use the NVECTOR_PARALLEL,
NVECTOR_PARHYP, NVECTOR_PETSC, or NVECTOR_TRILINOS implementations
underneath the NVECTOR_MPIPLUSX module since they already provide
MPI capabilities.


NVECTOR_MPIPLUSX structure
-------------------------------

The NVECTOR_MPIPLUSX implementation is a thin wrapper around the
NVECTOR_MPIMANYVECTOR. Accordingly, it adopts the same content
structure as defined in :numref:`NVectors.MPIManyVector.structure`.

The header file to include when using this module is
``nvector_mpiplusx.h``. The installed module library to link against is
``libsundials_nvecmpiplusx.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

.. note::

   If SUNDIALS is configured with MPI disabled, then the
   mpiplusx library will not be built.  Furthermore, any user codes
   that include ``nvector_mpiplusx.h`` *must* be compiled using an
   MPI-aware compiler.


NVECTOR_MPIPLUSX functions
-------------------------------

The NVECTOR_MPIPLUSX module adopts all vector operations listed
in :numref:`NVectors.Ops`, from the NVECTOR_MPIMANYVECTOR (see
:numref:`NVectors.MPIManyVector`) except for
:c:func:`N_VGetArrayPointer`, and :c:func:`N_VSetArrayPointer`;
the module provides its own implementation of these functions that
call the local vector implementations. Therefore, the NVECTOR_MPIPLUSX
module implements all of the operations listed in the referenced
sections except for :c:func:`N_VScaleAddMultiVectorArray()`,
and :c:func:`N_VLinearCombinationVectorArray()`. Accordingly, it's
compatibility with the SUNDIALS direct solvers and preconditioners
depends on the local vector implementation.

The module NVECTOR_MPIPLUSX provides the following additional
user-callable routines:

.. c:function:: N_Vector N_VMake_MPIPlusX(MPI_Comm comm, N_Vector *local_vector, SUNContext sunctx)

   This function creates a MPIPlusX vector from an exisiting local
   (i.e. on node) NVECTOR object, and a user-created MPI communicator.

   The input *comm* should be this user-created MPI communicator.
   This routine will internally call ``MPI_Comm_dup`` to create a
   copy of the input ``comm``, so the user-supplied ``comm`` argument
   need not be retained after the call to
   :c:func:`N_VMake_MPIPlusX`.

   This routine will copy the NVECTOR pointer to the input
   ``local_vector``, so the underlying local NVECTOR object
   should not be destroyed before the mpiplusx that contains it.

   Upon successful completion, the new MPIPlusX is returned;
   otherwise this routine returns ``NULL`` (e.g., if the input
   *local_vector* is ``NULL``).


.. c:function:: N_Vector N_VGetLocal_MPIPlusX(N_Vector v)

   This function returns the local vector underneath the MPIPlusX
   NVECTOR.


.. c:function:: realtype *N_VGetArrayPointer_MPIPlusX(N_Vector v)

   This function returns the data array pointer for the local vector.

   If the local vector does not support the :c:func:`N_VGetArrayPointer`
   operation, then ``NULL`` is returned.


.. c:function:: void N_VSetArrayPointer_MPIPlusX(realtype *v_data, N_Vector v)

   This function sets the data array pointer for the local vector if
   the local vector implements the :c:func:`N_VSetArrayPointer` operation.


The NVECTOR_MPIPLUSX module does not implement any fused or vector array
operations. Instead users should enable/disable fused operations on the
local vector.

**Notes**

* :c:func:`N_VMake_MPIPlusX` sets the field ``own_data = SUNFALSE`` and
  :c:func:`N_VDestroy_MPIPlusX()` will not call :c:func:`N_VDestroy()` on the
  local vector. In this a case, it is the user's responsibility to deallocate
  the local vector.

* To maximize efficiency, arithmetic vector operations in the
  NVECTOR_MPIPLUSX implementation that have more than one
  ``N_Vector`` argument do not check for consistent internal
  representation of these vectors. It is the user's responsibility to
  ensure that such routines are called with ``N_Vector`` arguments
  that were all created with the same subvector representations.
