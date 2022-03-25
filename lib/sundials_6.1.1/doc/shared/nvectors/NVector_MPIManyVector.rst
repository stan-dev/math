..
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

.. _NVectors.MPIManyVector:

The NVECTOR_MPIMANYVECTOR Module
================================

The NVECTOR_MPIMANYVECTOR module is designed to facilitate problems with an
inherent data partitioning for the solution vector, and when using
distributed-memory parallel architectures.  As such, this implementation
supports all use cases allowed by the MPI-unaware NVECTOR_MANYVECTOR
implementation, as well as partitioning data between nodes in a parallel
environment.  These data partitions are entirely user-defined, through
construction of distinct NVECTOR modules for each component, that are then
combined together to form the NVECTOR_MPIMANYVECTOR.  Three potential
use cases for this module include:

A. *Heterogenous computational architectures (single-node or multi-node)*:
   for data partitioning between different computing resources on a node,
   architecture-specific subvectors may be created for each partition.
   For example, a user could create one MPI-parallel component based on
   :ref:`NVECTOR_PARALLEL <NVectors.NVParallel>`, another GPU-accelerated
   component based on :ref:`NVECTOR_CUDA <NVectors.CUDA>`.

B. *Process-based multiphysics decompositions (multi-node)*:
   for computations that combine separate MPI-based simulations together,
   each subvector may reside on a different MPI communicator, and the
   MPIManyVector combines these via an MPI *intercommunicator* that
   connects these distinct simulations together.

C. *Structure of arrays (SOA) data layouts (single-node or multi-node)*:
   for problems that require
   separate subvectors for each solution component.  For example, in an
   incompressible Navier-Stokes simulation, separate subvectors may be
   used for velocities and pressure, which are combined together into a
   single MPIManyVector for the overall "solution".

The above use cases are neither exhaustive nor mutually exclusive, and
the NVECTOR_MANYVECTOR implementation should support arbitrary
combinations of these cases.

The NVECTOR_MPIMANYVECTOR implementation is designed to work with any
NVECTOR subvectors that implement the minimum "standard" set
of operations in :numref:`NVectors.Ops.Standard`, however significant
performance benefits may be obtained when subvectors additionally
implement the optional local reduction operations listed in
:numref:`NVectors.Ops.Local`.

Additionally, NVECTOR_MPIMANYVECTOR sets no limit on the number of
subvectors that may be attached (aside from the limitations of using
``sunindextype`` for indexing, and standard per-node memory
limitations).  However, while this ostensibly supports subvectors
with one entry each (i.e., one subvector for each solution entry), we
anticipate that this extreme situation will hinder performance due to
non-stride-one memory accesses and increased function call overhead.
We therefore recommend a relatively coarse partitioning of the
problem, although actual performance will likely be
problem-dependent.

As a final note, in the coming years we plan to introduce additional
algebraic solvers and time integration modules that will leverage the
problem partitioning enabled by NVECTOR_MPIMANYVECTOR.  However, even at
present we anticipate that users will be able to leverage such data
partitioning in their problem-defining ODE right-hand side function, DAE
or nonlinear solver residual function, preconditioners, or custom
:c:type:`SUNLinearSolver` or :c:type:`SUNNonlinearSolver` modules.


.. _NVectors.MPIManyVector.structure:

NVECTOR_MPIMANYVECTOR structure
-------------------------------

The NVECTOR_MPIMANYVECTOR implementation defines the *content* field
of ``N_Vector`` to be a structure containing the MPI communicator
(or ``MPI_COMM_NULL`` if running on a single-node), the number of
subvectors comprising the MPIManyVector, the global length of the
MPIManyVector (including all subvectors on all MPI ranks), a pointer to
the beginning of the array of subvectors, and a boolean flag
``own_data`` indicating ownership of the subvectors that populate
``subvec_array``.

.. code-block:: c

   struct _N_VectorContent_MPIManyVector {
     MPI_Comm      comm;            /* overall MPI communicator        */
     sunindextype  num_subvectors;  /* number of vectors attached      */
     sunindextype  global_length;   /* overall mpimanyvector length    */
     N_Vector*     subvec_array;    /* pointer to N_Vector array       */
     booleantype   own_data;        /* flag indicating data ownership  */
   };

The header file to include when using this module is
``nvector_mpimanyvector.h``. The installed module library to link against is
``libsundials_nvecmpimanyvector.lib`` where ``.lib`` is typically ``.so`` for
shared libraries and ``.a`` for static libraries.

.. note::

   If SUNDIALS is configured with MPI disabled, then the MPIManyVector
   library will not be built.  Furthermore, any user codes that include
   ``nvector_mpimanyvector.h`` *must* be compiled using an MPI-aware
   compiler (whether the specific user code utilizes MPI or not).  We
   note that the NVECTOR_MANYVECTOR implementation is designed for
   ManyVector use cases in an MPI-unaware environment.


NVECTOR_MPIMANYVECTOR functions
-------------------------------

The NVECTOR_MPIMANYVECTOR module implements all vector operations listed
in :numref:`NVectors.Ops`, except for :c:func:`N_VGetArrayPointer()`,
:c:func:`N_VSetArrayPointer()`, :c:func:`N_VScaleAddMultiVectorArray()`,
and :c:func:`N_VLinearCombinationVectorArray()`.  As such, this vector
cannot be used with the SUNDIALS direct solvers and preconditioners.
Instead, the NVECTOR_MPIMANYVECTOR module provides functions to access
subvectors, whose data may in turn be accessed according to their
NVECTOR implementations.

The names of vector operations are obtained from those in
:numref:`NVectors.Ops` by appending the suffix ``_MPIManyVector`` (e.g.
``N_VDestroy_MPIManyVector``).  The module NVECTOR_MPIMANYVECTOR provides
the following additional user-callable routines:

.. c:function:: N_Vector N_VNew_MPIManyVector(sunindextype num_subvectors, N_Vector *vec_array, SUNContext sunctx)

   This function creates a MPIManyVector from a set of existing
   NVECTOR objects, under the requirement that all MPI-aware
   subvectors use the same MPI communicator (this is checked
   internally).  If none of the subvectors are MPI-aware, then this
   may equivalently be used to describe data partitioning within a
   single node.  We note that this routine is designed to support use
   cases A and C above.

   This routine will copy all ``N_Vector`` pointers from the input
   ``vec_array``, so the user may modify/free that pointer array
   after calling this function.  However, this routine does *not*
   allocate any new subvectors, so the underlying NVECTOR objects
   themselves should not be destroyed before the MPIManyVector that
   contains them.

   Upon successful completion, the new MPIManyVector is returned;
   otherwise this routine returns ``NULL`` (e.g., if two MPI-aware
   subvectors use different MPI communicators).

   Users of the Fortran 2003 interface to this function will first need to use
   the generic ``N_Vector`` utility functions :c:func:`N_VNewVectorArray`, and
   :c:func:`N_VSetVecAtIndexVectorArray` to create the ``N_Vector*`` argument.  This is
   further explained in :numref:`SUNDIALS.Fortran.Differences.NVectorArrays`,
   and the functions are documented in :numref:`NVectors.Description.utilities`.


.. c:function:: N_Vector N_VMake_MPIManyVector(MPI_Comm comm, sunindextype num_subvectors, N_Vector *vec_array, SUNContext sunctx)

   This function creates a MPIManyVector from a set of existing NVECTOR
   objects, and a user-created MPI communicator that "connects" these
   subvectors.  Any MPI-aware subvectors may use different MPI
   communicators than the input *comm*.  We note that this routine
   is designed to support any combination of the use cases above.

   The input *comm* should be this user-created MPI communicator.
   This routine will internally call ``MPI_Comm_dup`` to create a
   copy of the input ``comm``, so the user-supplied ``comm`` argument
   need not be retained after the call to
   :c:func:`N_VMake_MPIManyVector`.

   If all subvectors are MPI-unaware, then the input *comm* argument
   should be ``MPI_COMM_NULL``, although in this case, it would be
   simpler to call :c:func:`N_VNew_MPIManyVector` instead, or to just
   use the NVECTOR_MANYVECTOR module.

   This routine will copy all ``N_Vector`` pointers from the input
   *vec_array*, so the user may modify/free that pointer array
   after calling this function.  However, this routine does *not*
   allocate any new subvectors, so the underlying NVECTOR objects
   themselves should not be destroyed before the MPIManyVector that
   contains them.

   Upon successful completion, the new MPIManyVector is returned;
   otherwise this routine returns ``NULL`` (e.g., if the input
   *vec_array* is ``NULL``).


.. c:function:: N_Vector N_VGetSubvector_MPIManyVector(N_Vector v, sunindextype vec_num)

   This function returns the *vec_num* subvector from the NVECTOR array.


.. c:function:: realtype *N_VGetSubvectorArrayPointer_MPIManyVector(N_Vector v, sunindextype vec_num)

   This function returns the data array pointer for the *vec_num*
   subvector from the NVECTOR array.

   If the input *vec_num* is invalid, or if the subvector does not
   support the ``N_VGetArrayPointer`` operation, then ``NULL`` is
   returned.


.. c:function:: int N_VSetSubvectorArrayPointer_MPIManyVector(realtype *v_data, N_Vector v, sunindextype vec_num)

   This function sets the data array pointer for the *vec_num*
   subvector from the NVECTOR array.

   If the input *vec_num* is invalid, or if the subvector does not
   support the ``N_VSetArrayPointer`` operation, then ``-1`` is
   returned; otherwise it returns ``0``.


.. c:function:: sunindextype N_VGetNumSubvectors_MPIManyVector(N_Vector v)

   This function returns the overall number of subvectors in the MPIManyVector object.


By default all fused and vector array operations are disabled in the
NVECTOR_MPIMANYVECTOR module, except for :c:func:`N_VWrmsNormVectorArray()`
and :c:func:`N_VWrmsNormMaskVectorArray()`, that are enabled by default.
The following additional user-callable routines are provided to enable or
disable fused and vector array operations for a specific vector. To
ensure consistency across vectors it is recommended to first create a
vector with :c:func:`N_VNew_MPIManyVector` or
:c:func:`N_VMake_MPIManyVector`, enable/disable the desired operations
for that vector with the functions below, and create any additional
vectors from that vector using :c:func:`N_VClone()`. This guarantees
that the new vectors will have the same operations enabled/disabled,
since cloned vectors inherit those configuration options from the
vector they are cloned from, while vectors created with
:c:func:`N_VNew_MPIManyVector` and :c:func:`N_VMake_MPIManyVector` will
have the default settings for the NVECTOR_MPIMANYVECTOR module.  We note
that these routines *do not* call the corresponding routines on
subvectors, so those should be set up as desired *before* attaching
them to the MPIManyVector in :c:func:`N_VNew_MPIManyVector` or
:c:func:`N_VMake_MPIManyVector`.

.. c:function:: int N_VEnableFusedOps_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) all fused and
   vector array operations in the MPIManyVector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableLinearCombination_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear
   combination fused operation in the MPIManyVector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleAddMulti_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale and
   add a vector to multiple vectors fused operation in the MPIManyVector vector. The
   return value is ``0`` for success and ``-1`` if the input vector or its
   ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableDotProdMulti_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the multiple
   dot products fused operation in the MPIManyVector vector. The return value is ``0``
   for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.

.. c:function:: int N_VEnableLinearSumVectorArray_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the linear sum
   operation for vector arrays in the MPIManyVector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableScaleVectorArray_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the scale
   operation for vector arrays in the MPIManyVector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableConstVectorArray_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the const
   operation for vector arrays in the MPIManyVector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormVectorArray_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the WRMS norm
   operation for vector arrays in the MPIManyVector vector. The return value is ``0`` for
   success and ``-1`` if the input vector or its ``ops`` structure are ``NULL``.

.. c:function:: int N_VEnableWrmsNormMaskVectorArray_MPIManyVector(N_Vector v, booleantype tf)

   This function enables (``SUNTRUE``) or disables (``SUNFALSE``) the masked WRMS
   norm operation for vector arrays in the MPIManyVector vector. The return value is
   ``0`` for success and ``-1`` if the input vector or its ``ops`` structure are
   ``NULL``.


**Notes**

* :c:func:`N_VNew_MPIManyVector` and :c:func:`N_VMake_MPIManyVector` set
  the field ``own_data = SUNFALSE``.
  :c:func:`N_VDestroy_MPIManyVector()` will not attempt to call
  :c:func:`N_VDestroy()` on any subvectors contained in the
  subvector array for any ``N_Vector`` with ``own_data`` set to
  ``SUNFALSE``. In such a case, it is the user's responsibility to
  deallocate the subvectors.

* To maximize efficiency, arithmetic vector operations in the
  NVECTOR_MPIMANYVECTOR implementation that have more than one
  ``N_Vector`` argument do not check for consistent internal
  representation of these vectors. It is the user's responsibility to
  ensure that such routines are called with ``N_Vector`` arguments
  that were all created with the same subvector representations.
