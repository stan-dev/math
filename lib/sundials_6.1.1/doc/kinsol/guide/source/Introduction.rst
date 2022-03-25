.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _KINSOL.Introduction:

************
Introduction
************

KINSOL is part of a software family called SUNDIALS: SUite of Nonlinear and DIfferential/ALgebraic equation
Solvers :cite:p:`HBGLSSW:05`. This suite consists of CVODE, ARKODE, KINSOL, and IDA,
and variants of these with sensitivity analysis capabilities.

KINSOL is a general-purpose nonlinear system solver based on Newton-Krylov solver technology. A fixed point
iteration is also included with the release of KINSOL v.2.8.0 and higher.

.. _KINSOL.Introduction.Historical:

Historical Background
=====================

The first nonlinear solver packages based on Newton-Krylov methods were written in Fortran. In particular, the NKSOL
package, written at LLNL, was the first Newton-Krylov solver package written for solution of systems arising in the
solution of partial differential equations :cite:p:`BrSa:90`. This Fortran code made use of Newton’s method to
solve the discrete nonlinear systems and applied a preconditioned Krylov linear solver for solution of the Jacobian
system at each nonlinear iteration. The key to the Newton-Krylov method was that the matrix-vector multiplies required
by the Krylov method could effectively be approximated by a finite difference of the nonlinear system-defining function,
avoiding a requirement for the formation of the actual Jacobian matrix. Significantly less memory was required for the
solver as a result.

In the late 1990’s, there was a push at LLNL to rewrite the nonlinear solver in C and port it to distributed
memory parallel machines. Both Newton and Krylov methods are easily implemented in parallel, and this effort gave rise
to the KINSOL package. KINSOL is similar to NKSOL in functionality, except that it provides for more options
in the choice of linear system methods and tolerances, and has a more modular design to provide flexibility for future
enhancements.

At present, KINSOL may utilize a variety of Krylov methods provided in SUNDIALS. These methods include the
GMRES (Generalized Minimal RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized Minimum
RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate Gradient Stabilized) :cite:p:`Van:92`, TFQMR
(Transpose-Free Quasi-Minimal Residual) :cite:p:`Fre:93`, and PCG (Preconditioned Conjugate
Gradient) :cite:p:`HeSt:52` linear iterative methods. As Krylov methods, these require little matrix storage
for solving the Newton equations as compared to direct methods. However, the algorithms allow for a user-supplied
preconditioner, and, for most problems, preconditioning is essential for an efficient solution. For very large
nonlinear algebraic systems, the Krylov methods are preferable over direct linear solver methods, and are often the only
feasible choice. Among the Krylov methods in SUNDIALS, we recommend GMRES as the best overall choice. However,
users are encouraged to compare all options, especially if encountering convergence failures with GMRES. Bi-CGStab and
TFQMR have an advantage in storage requirements, in that the number of workspace vectors they require is fixed, while
that number for GMRES depends on the desired Krylov subspace size. FGMRES has an advantage in that it is designed to
support preconditioners that vary between iterations (e.g., iterative methods). PCG exhibits rapid convergence and
minimal workspace vectors, but only works for symmetric linear systems.

For the sake of completeness in functionality, direct linear system solvers are included in KINSOL. These include
methods for both dense and banded linear systems, with Jacobians that are either user-supplied or generated internally
by difference quotients. KINSOL also includes interfaces to sparse direct solvers, including
KLU :cite:p:`DaPa:10,KLU_site` and the threaded sparse direct solver,
SuperLU_MT :cite:p:`Li:05,DGL:99,SuperLUMT_site`, among others (see Chapter :numref:`SUNLinSol` for further details).

In the process of translating NKSOL into C, the overall KINSOL organization has been changed considerably.
One key feature of the KINSOL organization is that a separate module devoted to vector operations was created.
This module facilitated extension to multiprosessor environments with minimal impact on the rest of the solver. The
vector module design is shared across the SUNDIALS suite. This :c:type:`N_Vector` module is written in terms of
abstract vector operations with the actual routines attached by a particular implementation (such as serial or parallel)
of ``N_Vector``. This abstraction allows writing the SUNDIALS solvers in a manner independent of the actual
``N_Vector`` implementation (which can be user-supplied), as well as allowing more than one ``N_Vector`` module linked
into an executable file. SUNDIALS (and thus KINSOL) is supplied with serial, MPI-parallel, OpenMP
and Pthreads thread-parallel ``N_Vector`` implementations, as well as multiple ``N_Vector``
implementations designed to leverage GPU architectures (see Chapter :numref:`NVectors` for
further details).

There are several motivations for choosing the C language for KINSOL. First, a general movement away from
Fortran and toward C in scientific computing was apparent. Second, the pointer, structure, and dynamic memory
allocation features in C are extremely useful in software of this complexity, with the great variety of method options
offered. Finally, we prefer C over C++ for KINSOL because of the wider availability of C
compilers, the potentially greater efficiency of C, and the greater ease of interfacing the solver to
applications written in Fortran.

.. _KINSOL.Introduction.Changes:

Changes from previous versions
==============================

Changes in v6.1.1
-----------------

Fixed exported ``SUNDIALSConfig.cmake``.

Changes in v6.1.0
-----------------

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the :c:func:`N_VSetKernelExecPolicy_Cuda` and
:c:func:`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

``SUNDIALS::<lib>`` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

:cmakeop:`CMAKE_C_STANDARD` is now set to 99 by default.

Fixed exported ``SUNDIALSConfig.cmake`` when profiling is enabled without Caliper.

Fixed ``sundials_export.h`` include in ``sundials_config.h``.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

Changes in v6.0.0
-----------------

**SUNContext**

SUNDIALS v6.0.0 introduces a new :c:type:`SUNContext` object on which all other
SUNDIALS objects depend. As such, the constructors for all SUNDIALS packages,
vectors, matrices, linear solvers, nonlinear solvers, and memory helpers have
been updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call :c:func:`SUNContext_Create` to create a context object
with before calling any other SUNDIALS library function, and then provide this
object to other SUNDIALS constructors. The context object has been introduced to
allow SUNDIALS to provide new features, such as the profiling/instrumentation
also introduced in this release, while maintaining thread-safety. See the
documentation section on the :c:type:`SUNContext` for more details.

A script ``upgrade-to-sundials-6-from-5.sh`` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a ``SUNCTX_PLACEHOLDER`` argument to all of
the calls to SUNDIALS constructors that now require a ``SUNContext`` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

.. code-block::

   > ./upgrade-to-sundials-6-from-5.sh <files to update>

**SUNProfiler**

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING`. A
built-in profiler will be used by default, but the `Caliper
<https://github.com/LLNL/Caliper>`_ library can also be used instead with the
CMake option :cmakeop:`ENABLE_CALIPER`. See the documentation section on
profiling for more details.  **WARNING**: Profiling will impact performance, and
should be enabled judiciously.

**SUNMemoryHelper**

The :c:type:`SUNMemoryHelper` functions :c:func:`SUNMemoryHelper_Alloc`,
:c:func:`SUNMemoryHelper_Dealloc`, and :c:func:`SUNMemoryHelper_Copy` have been
updated to accept an opaque handle as the last input. At a minimum, user-defined
:c:type:`SUNMemoryHelper` implementations will need to update these functions to
accept the additional argument. Typically, this handle is the execution stream
(e.g., a CUDA/HIP stream or SYCL queue) for the operation. The :ref:`CUDA
<SUNMemory.CUDA>`, :ref:`HIP <SUNMemory.HIP>`, and :ref:`SYCL <SUNMemory.SYCL>`
implementations have been updated accordingly. Additionally, the constructor
:c:func:`SUNMemoryHelper_Sycl` has been updated to remove the SYCL queue as an
input.

**NVector**

Two new optional vector operations, :c:func:`N_VDotProdMultiLocal` and
:c:func:`N_VDotProdMultiAllReduce`, have been added to support
low-synchronization methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the ``sundials``
namespace to the ``sundials::cuda``, ``sundials::hip``, and ``sundials::sycl``
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The ``Sundials`` namespace used by the Trilinos Tpetra NVector has been replaced
with the ``sundials::trilinos::nvector_tpetra`` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions ``N_VCloneVectorArray_*`` and ``N_VDestroyVectorArray_*`` have been
deprecated. The generic :c:func:`N_VCloneVectorArray` and
:c:func:`N_VDestroyVectorArray` functions should be used instead.

The previously deprecated constructor ``N_VMakeWithManagedAllocator_Cuda`` and
the function ``N_VSetCudaStream_Cuda`` have been removed and replaced with
:c:func:`N_VNewWithMemHelp_Cuda` and :c:func:`N_VSetKerrnelExecPolicy_Cuda`
respectively.

The previously deprecated macros ``PVEC_REAL_MPI_TYPE`` and
``PVEC_INTEGER_MPI_TYPE`` have been removed and replaced with
``MPI_SUNREALTYPE`` and ``MPI_SUNINDEXTYPE`` respectively.

**SUNLinearSolver**

The following previously deprecated functions have been removed:

+-----------------------------+------------------------------------------+
| Removed                     | Replacement                              |
+=============================+==========================================+
| ``SUNBandLinearSolver``     | :c:func:`SUNLinSol_Band`                 |
+-----------------------------+------------------------------------------+
| ``SUNDenseLinearSolver``    | :c:func:`SUNLinSol_Dense`                |
+-----------------------------+------------------------------------------+
| ``SUNKLU``                  | :c:func:`SUNLinSol_KLU`                  |
+-----------------------------+------------------------------------------+
| ``SUNKLUReInit``            | :c:func:`SUNLinSol_KLUReInit`            |
+-----------------------------+------------------------------------------+
| ``SUNKLUSetOrdering``       | :c:func:`SUNLinSol_KLUSetOrdering`       |
+-----------------------------+------------------------------------------+
| ``SUNLapackBand``           | :c:func:`SUNLinSol_LapackBand`           |
+-----------------------------+------------------------------------------+
| ``SUNLapackDense``          | :c:func:`SUNLinSol_LapackDense`          |
+-----------------------------+------------------------------------------+
| ``SUNPCG``                  | :c:func:`SUNLinSol_PCG`                  |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetPrecType``       | :c:func:`SUNLinSol_PCGSetPrecType`       |
+-----------------------------+------------------------------------------+
| ``SUNPCGSetMaxl``           | :c:func:`SUNLinSol_PCGSetMaxl`           |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGS``               | :c:func:`SUNLinSol_SPBCGS`               |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetPrecType``    | :c:func:`SUNLinSol_SPBCGSSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPBCGSSetMaxl``        | :c:func:`SUNLinSol_SPBCGSSetMaxl`        |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMR``               | :c:func:`SUNLinSol_SPFGMR`               |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetPrecType``    | :c:func:`SUNLinSol_SPFGMRSetPrecType`    |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetGSType``      | :c:func:`SUNLinSol_SPFGMRSetGSType`      |
+-----------------------------+------------------------------------------+
| ``SUNSPFGMRSetMaxRestarts`` | :c:func:`SUNLinSol_SPFGMRSetMaxRestarts` |
+-----------------------------+------------------------------------------+
| ``SUNSPGMR``                | :c:func:`SUNLinSol_SPGMR`                |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetPrecType``     | :c:func:`SUNLinSol_SPGMRSetPrecType`     |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetGSType``       | :c:func:`SUNLinSol_SPGMRSetGSType`       |
+-----------------------------+------------------------------------------+
| ``SUNSPGMRSetMaxRestarts``  | :c:func:`SUNLinSol_SPGMRSetMaxRestarts`  |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMR``              | :c:func:`SUNLinSol_SPTFQMR`              |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetPrecType``   | :c:func:`SUNLinSol_SPTFQMRSetPrecType`   |
+-----------------------------+------------------------------------------+
| ``SUNSPTFQMRSetMaxl``       | :c:func:`SUNLinSol_SPTFQMRSetMaxl`       |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMT``            | :c:func:`SUNLinSol_SuperLUMT`            |
+-----------------------------+------------------------------------------+
| ``SUNSuperLUMTSetOrdering`` | :c:func:`SUNLinSol_SuperLUMTSetOrdering` |
+-----------------------------+------------------------------------------+

**KINSOL**

New orthogonalization methods were added for use within the KINSOL Anderson
acceleration routine. See :numref:`Anderson_QR` and :c:func:`KINSetOrthAA`
for more details.

The KINSOL Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
and the F2003 example programs for more details using the SUNDIALS Fortran 2003
module interfaces.

**Deprecations**

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and typedefs are now deprecated:

+------------------------------+-------------------------------------+
| Deprecated Name              | New Name                            |
+==============================+=====================================+
| ``realtype``                 | ``sunrealtype``                     |
+------------------------------+-------------------------------------+
| ``booleantype``              | ``sunbooleantype``                  |
+------------------------------+-------------------------------------+
| ``RCONST``                   | ``SUN_RCONST``                      |
+------------------------------+-------------------------------------+
| ``BIG_REAL``                 | ``SUN_BIG_REAL``                    |
+------------------------------+-------------------------------------+
| ``SMALL_REAL``               | ``SUN_SMALL_REAL``                  |
+------------------------------+-------------------------------------+
| ``UNIT_ROUNDOFF``            | ``SUN_UNIT_ROUNDOFF``               |
+------------------------------+-------------------------------------+
| ``PREC_NONE``                | ``SUN_PREC_NONE``                   |
+------------------------------+-------------------------------------+
| ``PREC_LEFT``                | ``SUN_PREC_LEFT``                   |
+------------------------------+-------------------------------------+
| ``PREC_RIGHT``               | ``SUN_PREC_RIGHT``                  |
+------------------------------+-------------------------------------+
| ``PREC_BOTH``                | ``SUN_PREC_BOTH``                   |
+------------------------------+-------------------------------------+
| ``MODIFIED_GS``              | ``SUN_MODIFIED_GS``                 |
+------------------------------+-------------------------------------+
| ``CLASSICAL_GS``             | ``SUN_CLASSICAL_GS``                |
+------------------------------+-------------------------------------+
| ``ATimesFn``                 | ``SUNATimesFn``                     |
+------------------------------+-------------------------------------+
| ``PSetupFn``                 | ``SUNPSetupFn``                     |
+------------------------------+-------------------------------------+
| ``PSolveFn``                 | ``SUNPSolveFn``                     |
+------------------------------+-------------------------------------+
| ``DlsMat``                   | ``SUNDlsMat``                       |
+------------------------------+-------------------------------------+
| ``DENSE_COL``                | ``SUNDLS_DENSE_COL``                |
+------------------------------+-------------------------------------+
| ``DENSE_ELEM``               | ``SUNDLS_DENSE_ELEM``               |
+------------------------------+-------------------------------------+
| ``BAND_COL``                 | ``SUNDLS_BAND_COL``                 |
+------------------------------+-------------------------------------+
| ``BAND_COL_ELEM``            | ``SUNDLS_BAND_COL_ELEM``            |
+------------------------------+-------------------------------------+
| ``BAND_ELEM``                | ``SUNDLS_BAND_ELEM``                |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
| ``KINDlsSetLinearSolver``       | ``KINSetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``KINDlsSetJacFn``              | ``KINSetJacFn``                |
+---------------------------------+--------------------------------+
| ``KINDlsGetWorkSpace``          | ``KINGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``KINDlsGetNumJacEvals``        | ``KINGetNumJacEvals``          |
+---------------------------------+--------------------------------+
| ``KINDlsGetNumFuncEvals``       | ``KINGetNumLinFuncEvals``      |
+---------------------------------+--------------------------------+
| ``KINDlsGetLastFlag``           | ``KINGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``KINDlsGetReturnFlagName``     | ``KINGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``KINSpilsSetLinearSolver``     | ``KINSetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``KINSpilsSetPreconditioner``   | ``KINSetPreconditioner``       |
+---------------------------------+--------------------------------+
| ``KINSpilsSetJacTimesVecFn``    | ``KINSetJacTimesVecFn``        |
+---------------------------------+--------------------------------+
| ``KINSpilsGetWorkSpace``        | ``KINGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumPrecEvals``     | ``KINGetNumPrecEvals``         |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumPrecSolves``    | ``KINGetNumPrecSolves``        |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumLinIters``      | ``KINGetNumLinIters``          |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumConvFails``     | ``KINGetNumLinConvFails``      |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumJtimesEvals``   | ``KINGetNumJtimesEvals``       |
+---------------------------------+--------------------------------+
| ``KINSpilsGetNumFuncEvals``     | ``KINGetNumLinFuncEvals``      |
+---------------------------------+--------------------------------+
| ``KINSpilsGetLastFlag``         | ``KINGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``KINSpilsGetReturnFlagName``   | ``KINGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``DenseGETRF``                  | ``SUNDlsMat_DenseGETRF``       |
+---------------------------------+--------------------------------+
| ``DenseGETRS``                  | ``SUNDlsMat_DenseGETRS``       |
+---------------------------------+--------------------------------+
| ``denseGETRF``                  | ``SUNDlsMat_denseGETRF``       |
+---------------------------------+--------------------------------+
| ``denseGETRS``                  | ``SUNDlsMat_denseGETRS``       |
+---------------------------------+--------------------------------+
| ``DensePOTRF``                  | ``SUNDlsMat_DensePOTRF``       |
+---------------------------------+--------------------------------+
| ``DensePOTRS``                  | ``SUNDlsMat_DensePOTRS``       |
+---------------------------------+--------------------------------+
| ``densePOTRF``                  | ``SUNDlsMat_densePOTRF``       |
+---------------------------------+--------------------------------+
| ``densePOTRS``                  | ``SUNDlsMat_densePOTRS``       |
+---------------------------------+--------------------------------+
| ``DenseGEQRF``                  | ``SUNDlsMat_DenseGEQRF``       |
+---------------------------------+--------------------------------+
| ``DenseORMQR``                  | ``SUNDlsMat_DenseORMQR``       |
+---------------------------------+--------------------------------+
| ``denseGEQRF``                  | ``SUNDlsMat_denseGEQRF``       |
+---------------------------------+--------------------------------+
| ``denseORMQR``                  | ``SUNDlsMat_denseORMQR``       |
+---------------------------------+--------------------------------+
| ``DenseCopy``                   | ``SUNDlsMat_DenseCopy``        |
+---------------------------------+--------------------------------+
| ``denseCopy``                   | ``SUNDlsMat_denseCopy``        |
+---------------------------------+--------------------------------+
| ``DenseScale``                  | ``SUNDlsMat_DenseScale``       |
+---------------------------------+--------------------------------+
| ``denseScale``                  | ``SUNDlsMat_denseScale``       |
+---------------------------------+--------------------------------+
| ``denseAddIdentity``            | ``SUNDlsMat_denseAddIdentity`` |
+---------------------------------+--------------------------------+
| ``DenseMatvec``                 | ``SUNDlsMat_DenseMatvec``      |
+---------------------------------+--------------------------------+
| ``denseMatvec``                 | ``SUNDlsMat_denseMatvec``      |
+---------------------------------+--------------------------------+
| ``BandGBTRF``                   | ``SUNDlsMat_BandGBTRF``        |
+---------------------------------+--------------------------------+
| ``bandGBTRF``                   | ``SUNDlsMat_bandGBTRF``        |
+---------------------------------+--------------------------------+
| ``BandGBTRS``                   | ``SUNDlsMat_BandGBTRS``        |
+---------------------------------+--------------------------------+
| ``bandGBTRS``                   | ``SUNDlsMat_bandGBTRS``        |
+---------------------------------+--------------------------------+
| ``BandCopy``                    | ``SUNDlsMat_BandCopy``         |
+---------------------------------+--------------------------------+
| ``bandCopy``                    | ``SUNDlsMat_bandCopy``         |
+---------------------------------+--------------------------------+
| ``BandScale``                   | ``SUNDlsMat_BandScale``        |
+---------------------------------+--------------------------------+
| ``bandScale``                   | ``SUNDlsMat_bandScale``        |
+---------------------------------+--------------------------------+
| ``bandAddIdentity``             | ``SUNDlsMat_bandAddIdentity``  |
+---------------------------------+--------------------------------+
| ``BandMatvec``                  | ``SUNDlsMat_BandMatvec``       |
+---------------------------------+--------------------------------+
| ``bandMatvec``                  | ``SUNDlsMat_bandMatvec``       |
+---------------------------------+--------------------------------+
| ``ModifiedGS``                  | ``SUNModifiedGS``              |
+---------------------------------+--------------------------------+
| ``ClassicalGS``                 | ``SUNClassicalGS``             |
+---------------------------------+--------------------------------+
| ``QRfact``                      | ``SUNQRFact``                  |
+---------------------------------+--------------------------------+
| ``QRsol``                       | ``SUNQRsol``                   |
+---------------------------------+--------------------------------+
| ``DlsMat_NewDenseMat``          | ``SUNDlsMat_NewDenseMat``      |
+---------------------------------+--------------------------------+
| ``DlsMat_NewBandMat``           | ``SUNDlsMat_NewBandMat``       |
+---------------------------------+--------------------------------+
| ``DestroyMat``                  | ``SUNDlsMat_DestroyMat``       |
+---------------------------------+--------------------------------+
| ``NewIntArray``                 | ``SUNDlsMat_NewIntArray``      |
+---------------------------------+--------------------------------+
| ``NewIndexArray``               | ``SUNDlsMat_NewIndexArray``    |
+---------------------------------+--------------------------------+
| ``NewRealArray``                | ``SUNDlsMat_NewRealArray``     |
+---------------------------------+--------------------------------+
| ``DestroyArray``                | ``SUNDlsMat_DestroyArray``     |
+---------------------------------+--------------------------------+
| ``AddIdentity``                 | ``SUNDlsMat_AddIdentity``      |
+---------------------------------+--------------------------------+
| ``SetToZero``                   | ``SUNDlsMat_SetToZero``        |
+---------------------------------+--------------------------------+
| ``PrintMat``                    | ``SUNDlsMat_PrintMat``         |
+---------------------------------+--------------------------------+
| ``newDenseMat``                 | ``SUNDlsMat_newDenseMat``      |
+---------------------------------+--------------------------------+
| ``newBandMat``                  | ``SUNDlsMat_newBandMat``       |
+---------------------------------+--------------------------------+
| ``destroyMat``                  | ``SUNDlsMat_destroyMat``       |
+---------------------------------+--------------------------------+
| ``newIntArray``                 | ``SUNDlsMat_newIntArray``      |
+---------------------------------+--------------------------------+
| ``newIndexArray``               | ``SUNDlsMat_newIndexArray``    |
+---------------------------------+--------------------------------+
| ``newRealArray``                | ``SUNDlsMat_newRealArray``     |
+---------------------------------+--------------------------------+
| ``destroyArray``                | ``SUNDlsMat_destroyArray``     |
+---------------------------------+--------------------------------+

In addition, the entire ``sundials_lapack.h`` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

Changes in v5.8.0
-----------------

The RAJA ``N_Vector`` implementation has been updated to support the SYCL backend in addition to the CUDA and HIP
backend. Users can choose the backend when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable.
This module remains experimental and is subject to change from version to version.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to interface with the Intel oneAPI Math Kernel
Library (oneMKL). Both the matrix and the linear solver support general dense linear systems as well as block diagonal
linear systems. See :numref:`SUNLinSol.OneMklDense`  for more details. This module is
experimental and is subject to change from version to version.

Added a new *optional* function to the SUNLinearSolver API, ``SUNLinSolSetZeroGuess``, to indicate that the next call to
``SUNlinSolSolve`` will be made with a zero initial guess. SUNLinearSolver implementations that do not use the
``SUNLinSolNewEmpty`` constructor will, at a minimum, need set the ``setzeroguess`` function pointer in the linear
solver ``ops`` structure to ``NULL``. The SUNDIALS iterative linear solver implementations have been updated to leverage
this new set function to remove one dot product per solve.

New KINSOL options have been added to apply a constant damping in the fixed point and Picard iterations (see
``KINSetDamping``), to delay the start of Anderson acceleration with the fixed point and Picard iterations (see
``KINSetDelayAA``), and to return the newest solution with the fixed point iteration (see ``KINSetReturnNewest``).

The installed SUNDIALSConfig.cmake file now supports the ``COMPONENTS`` option to ``find_package``. The exported targets
no longer have ``IMPORTED_GLOBAL`` set.

A bug was fixed in ``SUNMatCopyOps`` where the matrix-vector product setup function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero initial guess and a solution scaling
vector are provided. This fix only impacts codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS
packages utilize a zero initial guess.

A bug was fixed in the Picard iteration where the value of ``KINSetMaxSetupCalls`` would be ignored.

Changes in v5.7.0
-----------------

A new ``N_Vector`` implementation based on the SYCL abstraction layer has been added targeting Intel GPUs. At
present the only SYCL compiler supported is the DPC++ (Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details. This module is considered experimental and is subject to major
changes even in minor releases.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to interface with the MAGMA linear algebra library.
Both the matrix and the linear solver support general dense linear systems as well as block diagonal linear systems, and
both are targeted at GPUs (AMD or NVIDIA). See :numref:`SUNLinSol.MagmaDense` for more
details.

Changes in v5.6.1
-----------------

Fixed a bug in the SUNDIALS CMake which caused an error if the CMAKE_CXX_STANDARD and SUNDIALS_RAJA_BACKENDS
options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

Changes in v5.6.0
-----------------

A new ``N_Vector`` implementation based on the AMD ROCm HIP platform has been added. This vector can target NVIDIA or
AMD GPUs. See :numref:`NVectors.HIP` for more details. This module is considered experimental and is subject
to change from version to version.

The RAJA ``N_Vector`` implementation has been updated to support the HIP backend in addition to the CUDA backend. Users
can choose the backend when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, ``N_VGetDeviceArrayPointer``, was added to the N_Vector API. This operation is useful for
N_Vectors that utilize dual memory spaces, e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations no longer require the SUNDIALS CUDA
N_Vector. Instead, they require that the vector utilized provides the ``N_VGetDeviceArrayPointer`` operation, and that
the pointer returned by ``N_VGetDeviceArrayPointer`` is a valid CUDA device pointer.

Changes in v5.5.0
-----------------

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required. Users will likely see deprecation
warnings, but otherwise the changes should be fully backwards compatible for almost all users. SUNDIALS now
exports CMake targets and installs a SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.

Changes in v5.4.0
-----------------

A new API, ``SUNMemoryHelper``, was added to support **GPU users** who have complex memory management needs such as
using memory pools. This is paired with new constructors for the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules that accept a
``SUNMemoryHelper`` object. Refer to :numref:`SUNDIALS.GPU.Model`, :numref:`NVectors.CUDA`, :numref:`NVectors.RAJA`, and :numref:`SUNMemory` for more information.

The ``NVECTOR_RAJA`` module has been updated to mirror the ``NVECTOR_CUDA`` module. Notably, the update adds managed
memory support to the ``NVECTOR_RAJA`` module. Users of the module will need to update any calls to the ``N_VMake_Raja``
function because that signature was changed. This module remains experimental and is subject to change from version to
version.

The ``NVECTOR_TRILINOS`` module has been updated to work with Trilinos 12.18+. This update changes the local ordinal
type to always be an ``int``.

Added support for CUDA v11.

Changes in v5.3.0
-----------------

Fixed a bug in the iterative linear solver modules where an error is not returned if the Atimes function is ``NULL`` or,
if preconditioning is enabled, the PSolve function is ``NULL``.

Added the ability to control the CUDA kernel launch parameters for the ``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE``
modules. These modules remain experimental and are subject to change from version to version. In addition, the
``NVECTOR_CUDA`` kernels were rewritten to be more flexible. Most users should see equivalent performance or some
improvement, but a select few may observe minor performance degradation with the default settings. Users are encouraged
to contact the SUNDIALS team about any perfomance changes that they notice.

Added new capabilities for monitoring the solve phase in the ``SUNNONLINSOL_NEWTON`` and ``SUNNONLINSOL_FIXEDPOINT``
modules, and the SUNDIALS iterative linear solver modules. SUNDIALS must be built with the CMake option
``SUNDIALS_BUILD_WITH_MONITORING`` to use these capabilties.

Added the optional function ``KINSetJacTimesVecSysFn`` to specify an alternative system function for computing
Jacobian-vector products with the internal difference quotient approximation.

Changes in v5.2.0
-----------------

Fixed a build system bug related to the Fortran 2003 interfaces when using the IBM XL compiler. When building the
Fortran 2003 interfaces with an XL compiler it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``,
``xlf2003``, or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from dllimport/dllexport attributes missing on some
SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, ``SUNMATRIX_CUSPARSE``, that interfaces to the sparse matrix implementation
from the NVIDIA cuSPARSE library. In addition, the ``SUNLINSOL_CUSOLVER_BATCHQR`` linear solver has been updated to use
this matrix, therefore, users of this module will need to update their code. These modules are still considered to be
experimental, thus they are subject to breaking changes even in minor releases.

Changes in v5.1.0
-----------------

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake variables ``PETSC_INCLUDES`` and
``PETSC_LIBRARIES`` instead of ``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify the CUDA architecture to compile for.

Added two utility functions, ``SUNDIALSFileOpen`` and ``SUNDIALSFileClose`` for creating/destroying file pointers that
are useful when using the Fortran 2003 interfaces.

Added support for constant damping when using Anderson acceleration. See :numref:`KINSOL.Mathematics` and the
description of the ``KINSetDampingAA`` function for more details.

Changes in v5.0.0
-----------------

Build system changes
^^^^^^^^^^^^^^^^^^^^

-  Increased the minimum required CMake version to 3.5 for most SUNDIALS configurations, and 3.10 when CUDA or
   OpenMP with device offloading are enabled.

-  The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been removed to simplify builds as
   SUNDIALS packages do not use BLAS directly. For third party libraries that require linking to BLAS, the path to
   the BLAS library should be included in the ``_LIBRARIES`` variable for the third party library *e.g.*,
   ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

-  Fixed a bug in the build system that prevented the ``NVECTOR_PTHREADS`` module from being built.

NVECTOR module changes
^^^^^^^^^^^^^^^^^^^^^^

-  Two new functions were added to aid in creating custom ``N_Vector`` objects. The constructor ``N_VNewEmpty``
   allocates an “empty” generic ``N_Vector`` with the object’s content pointer and the function pointers in the
   operations structure initialized to ``NULL``. When used in the constructor for custom objects this function will ease
   the introduction of any new optional operations to the ``N_Vector`` API by ensuring only required operations need to
   be set. Additionally, the function ``N_VCopyOps(w, v)`` has been added to copy the operation function pointers
   between vector objects. When used in clone routines for custom vector objects these functions also will ease the
   introduction of any new optional operations to the ``N_Vector`` API by ensuring all operations are copied when
   cloning objects. See :numref:`NVectors.Description.utilities` for more details.

-  Two new ``N_Vector`` implementations, ``NVECTOR_MANYVECTOR`` and ``NVECTOR_MPIMANYVECTOR``, have been created to support
   flexible partitioning of solution data among different processing elements (e.g., CPU + GPU) or for multi-physics
   problems that couple distinct MPI-based simulations together. This implementation is accompanied by additions to user
   documentation and SUNDIALS examples. See :numref:`NVectors.ManyVector` and :numref:`NVectors.MPIManyVector` for more details.

-  One new required vector operation and ten new optional vector operations have been added to the ``N_Vector`` API.
   The new required operation, ``N_VGetLength``, returns the global length of an ``N_Vector``. The optional operations
   have been added to support the new ``NVECTOR_MPIMANYVECTOR`` implementation. The operation ``N_VGetCommunicator`` must
   be implemented by subvectors that are combined to create an ``NVECTOR_MPIMANYVECTOR``, but is not used outside of this
   context. The remaining nine operations are optional local reduction operations intended to eliminate unnecessary
   latency when performing vector reduction operations (norms, etc.) on distributed memory systems. The optional local
   reduction vector operations are ``N_VDotProdLocal``, ``N_VMaxNormLocal``, ``N_VMinLocal``, ``N_VL1NormLocal``,
   ``N_VWSqrSumLocal``, ``N_VWSqrSumMaskLocal``, ``N_VInvTestLocal``, ``N_VConstrMaskLocal``, and
   ``N_VMinQuotientLocal``. If an ``N_Vector`` implementation defines any of the local operations as ``NULL``, then the
   ``NVECTOR_MPIMANYVECTOR`` will call standard ``N_Vector`` operations to complete the computation. See
   :numref:`NVectors.Ops.Local` for more details.

-  An additional ``N_Vector`` implementation, ``NVECTOR_MPIPLUSX``, has been created to support the MPI+X paradigm where
   X is a type of on-node parallelism (*e.g.*, OpenMP, CUDA). The implementation is accompanied by additions to user
   documentation and SUNDIALS examples. See :numref:`NVectors.MPIPlusX` for more details.

-  The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA``
   implementations respectively. Accordingly, the ``nvector_mpicuda.h``, ``nvector_mpiraja.h``,
   ``libsundials_nvecmpicuda.lib``, and ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should use
   the ``NVECTOR_MPIPLUSX`` module coupled in conjunction with the ``NVECTOR_CUDA`` or ``NVECTOR_RAJA`` modules to replace the
   functionality. The necessary changes are minimal and should require few code modifications. See the programs in
   ``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to use the ``NVECTOR_MPIPLUSX`` module with
   the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules respectively.

-  Fixed a memory leak in the ``NVECTOR_PETSC`` module clone function.

-  Made performance improvements to the ``NVECTOR_CUDA`` module. Users who utilize a non-default stream should no longer
   see default stream synchronizations after memory transfers.

-  Added a new constructor to the ``NVECTOR_CUDA`` module that allows a user to provide custom allocate and free functions
   for the vector data array and internal reduction buffer. See
   :numref:`NVectors.CUDA.Functions` for more details.

-  Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See Chapter :numref:`NVectors` for more
   details on how to use the interfaces.

-  Added three new ``N_Vector`` utility functions, ``FN_VGetVecAtIndexVectorArray``, ``FN_VSetVecAtIndexVectorArray``,
   and ``FN_VNewVectorArray``, for working with ``N_Vector`` arrays when using the Fortran
   2003 interfaces. See :numref:`NVectors.Description.utilities` for more details.

SUNMatrix module changes
^^^^^^^^^^^^^^^^^^^^^^^^

-  Two new functions were added to aid in creating custom ``SUNMatrix`` objects. The constructor ``SUNMatNewEmpty``
   allocates an “empty” generic ``SUNMatrix`` with the object’s content pointer and the function pointers in the
   operations structure initialized to ``NULL``. When used in the constructor for custom objects this function will ease
   the introduction of any new optional operations to the ``SUNMatrix`` API by ensuring only required operations need
   to be set. Additionally, the function ``SUNMatCopyOps(A, B)`` has been added to copy the operation function pointers
   between matrix objects. When used in clone routines for custom matrix objects these functions also will ease the
   introduction of any new optional operations to the ``SUNMatrix`` API by ensuring all operations are copied when
   cloning objects. See :numref:`SUNMatrix.Description` for more details.

-  A new operation, ``SUNMatMatvecSetup``, was added to the ``SUNMatrix`` API to perform any setup necessary for
   computing a matrix-vector product. This operation is useful for ``SUNMatrix`` implementations which need to prepare
   the matrix itself, or communication structures before performing the matrix-vector product. Users who have
   implemented custom ``SUNMatrix`` modules will need to at least update their code to set the corresponding ``ops``
   structure member, ``matvecsetup``, to ``NULL``. See :numref:`SUNMatrix.Ops` for
   more details.

-  The generic ``SUNMatrix`` API now defines error codes to be returned by ``SUNMatrix`` operations. Operations
   which return an integer flag indiciating success/failure may return different values than previously. See
   :numref:`SUNMatrix.Ops.errorCodes` for more details.

-  A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to facilitate the use of the SuperLU_DIST
   library with SUNDIALS. See :numref:`SUNMatrix.SLUNRloc` for more details.

-  Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See Chapter :numref:`SUNMatrix` for
   more details on how to use the interfaces.

SUNLinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  A new function was added to aid in creating custom ``SUNLinearSolver`` objects. The constructor ``SUNLinSolNewEmpty``
   allocates an “empty” generic ``SUNLinearSolver`` with the object’s content pointer and the function pointers in the
   operations structure initialized to ``NULL``. When used in the constructor for custom objects this function will ease
   the introduction of any new optional operations to the ``SUNLinearSolver`` API by ensuring only required operations need
   to be set. See :numref:`SUNLinSol.API.Custom` for more details.

-  The return type of the ``SUNLinearSolver`` API function ``SUNLinSolLastFlag`` has changed from ``long int`` to
   ``sunindextype`` to be consistent with the type used to store row indices in dense and banded linear solver modules.

-  Added a new optional operation to the ``SUNLinearSolver`` API, ``SUNLinSolGetID``, that returns a ``SUNLinearSolver_ID``
   for identifying the linear solver module.

-  The ``SUNLinearSolver`` API has been updated to make the initialize and setup functions optional.

-  A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to facilitate the use of the SuperLU_DIST
   library with SUNDIALS. See :numref:`SUNLinSol.SuperLUDIST` for more details.

-  Added a new ``SUNLinearSolver`` implementation, ``SUNLinearSolver_cuSolverSp_batchQR``, which leverages the NVIDIA
   cuSOLVER sparse batched QR method for efficiently solving block diagonal linear systems on NVIDIA GPUs. See :numref:`SUNLinSol.cuSolverSp` for more details.

-  Added three new accessor functions to the ``SUNLINSOL_KLU`` module, ``SUNLinSol_KLUGetSymbolic``,
   ``SUNLinSol_KLUGetNumeric``, and ``SUNLinSol_KLUGetCommon``, to provide user access to the underlying KLU solver
   structures. See :numref:`SUNLinSol.KLU.Usage` for more details.

-  Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules. See Chapter :numref:`SUNLinSol` for
   more details on how to use the interfaces.

KINSOL changes
^^^^^^^^^^^^^^

-  Fixed a bug in the KINSOL linear solver interface where the auxiliary scalar ``sJpnorm`` was not computed when
   necessary with the Picard iteration and the auxiliary scalar ``sFdotJp`` was unnecessarily computed in some cases.

-  The KINLS interface has been updated to only zero the Jacobian matrix before calling a user-supplied Jacobian
   evaluation function when the attached linear solver has type ``SUNLINEARSOLVER_DIRECT``.

-  Added a Fortran 2003 interface to KINSOL. See :numref:`SUNDIALS.Fortran` for more details.

Changes in v4.1.0
-----------------

An additional ``N_Vector`` implementation was added for the TPetra vector from the Trilinos library to
facilitate interoperability between SUNDIALS and Trilinos. This implementation is accompanied by additions
to user documentation and SUNDIALS examples.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option ``EXAMPLES_ENABLE_CUDA`` enables all examples
that use CUDA including the RAJA examples with a CUDA back end (if the RAJA ``N_Vector`` is enabled).

The implementation header file ``kin_impl.h`` is no longer installed. This means users who are directly manipulating the
``KINMem`` structure will need to update their code to use KINSOL’s public API.

Python is no longer required to run ``make test`` and ``make test_install``.

Changes in v4.0.2
-----------------

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source file. The symbols are now included in
the KINSOL library, ``libsundials_kinsol``.

Changes in v4.0.1
-----------------

No changes were made in this release.

Changes in v4.0.0
-----------------

KINSOL’s previous direct and iterative linear solver interfaces, KINDls and KINSpils, have been merged
into a single unified linear solver interface, KINLs, to support any valid ``SUNLinearSolver`` module. This includes
the “DIRECT” and “ITERATIVE” types as well as the new “MATRIX_ITERATIVE” type. Details regarding how KINLs
utilizes linear solvers of each type as well as discussion regarding intended use cases for user-supplied
``SUNLinearSolver`` implementations are included in Chapter :numref:`SUNLinSol`. All KINSOL example
programs and the standalone linear solver examples have been updated to use the unified linear solver interface.

The unified interface for the new KINLs module is very similar to the previous KINDls and KINSpils
interfaces. To minimize challenges in user migration to the new names, the previous C and Fortran routine names
may still be used; these will be deprecated in future releases, so we recommend that users migrate to the new names
soon. Additionally, we note that Fortran users, however, may need to enlarge their ``iout`` array of optional integer
outputs, and update the indices that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver`` implementations have been updated to
follow the naming convention ``SUNLinSol_*`` where ``*`` is the name of the linear solver. The new names are
``SUNLinSol_Band``, ``SUNLinSol_Dense``, ``SUNLinSol_KLU``, ``SUNLinSol_LapackBand``, ``SUNLinSol_LapackDense``,
``SUNLinSol_PCG``, ``SUNLinSol_SPBCGS``, ``SUNLinSol_SPFGMR``, ``SUNLinSol_SPGMR``, ``SUNLinSol_SPTFQMR``, and
``SUNLinSol_SuperLUMT``. Solver-specific “set” routine names have been similarly standardized. To minimize challenges in
user migration to the new names, the previous routine names may still be used; these will be deprecated in future
releases, so we recommend that users migrate to the new names soon. All KINSOL example programs and the standalone
linear solver examples have been updated to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the storage upper bandwidth argument.

Three fused vector operations and seven vector array operations have been added to the ``N_Vector`` API.
These *optional* operations are disabled by default and may be activated by calling vector specific
routines after creating an ``N_Vector`` (see Chapter :numref:`NVectors` for more details). The new
operations are intended to increase data reuse in vector operations, reduce parallel communication on
distributed memory systems, and lower the number of kernel launches on systems with accelerators. The
fused operations are ``N_VLinearCombination``, ``N_VScaleAddMulti``, and ``N_VDotProdMulti`` and the
vector array operations are ``N_VLinearCombinationVectorArray``, ``N_VScaleVectorArray``,
``N_VConstVectorArray``, ``N_VWrmsNormVectorArray``, ``N_VWrmsNormMaskVectorArray``,
``N_VScaleAddMultiVectorArray``, and ``N_VLinearCombinationVectorArray``. If an ``N_Vector``
implementation defines any of these operations as ``NULL``, then standard ``N_Vector`` operations will
automatically be called as necessary to complete the computation. Multiple updates to ``NVECTOR_CUDA``
were made:

-  Changed ``N_VGetLength_Cuda`` to return the global vector length instead of the local vector length.

-  Added ``N_VGetLocalLength_Cuda`` to return the local vector length.

-  Added ``N_VGetMPIComm_Cuda`` to return the MPI communicator used.

-  Removed the accessor functions in the namespace suncudavec.

-  Changed the ``N_VMake_Cuda`` function to take a host data pointer and a device data pointer instead of an
   ``N_VectorContent_Cuda`` object.

-  Added the ability to set the ``cudaStream_t`` used for execution of the ``NVECTOR_CUDA`` kernels. See the function
   ``N_VSetCudaStreams_Cuda``.

-  Added ``N_VNewManaged_Cuda``, ``N_VMakeManaged_Cuda``, and ``N_VIsManagedMemory_Cuda`` functions to accommodate using
   managed memory with the ``NVECTOR_CUDA``.

Multiple changes to ``NVECTOR_RAJA`` were made:

-  Changed ``N_VGetLength_Raja`` to return the global vector length instead of the local vector length.

-  Added ``N_VGetLocalLength_Raja`` to return the local vector length.

-  Added ``N_VGetMPIComm_Raja`` to return the MPI communicator used.

-  Removed the accessor functions in the namespace suncudavec.

A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device offloading has been added, ``NVECTOR_OPENMPDEV``. See
:numref:`NVectors.OpenMPDEV` for more details.

Changes in v3.2.1
-----------------

The changes in this minor release include the following:

-  Fixed a bug in the CUDA ``N_Vector`` where the ``N_VInvTest`` operation could write beyond the allocated
   vector data.

-  Fixed library installation path for multiarch systems. This fix changes the default library
   installation path to ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from ``CMAKE_INSTALL_PREFIX/lib``.
   ``CMAKE_INSTALL_LIBDIR`` is automatically set, but is available as a CMake option that can modified.

Changes in v3.2.0
-----------------

Fixed a problem with setting ``sunindextype`` which would occur with some compilers (e.g. armclang) that
did not define ``__STDC_VERSION__``. Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than
one MPI rank when using a GPU system. The vectors assume one GPU device per MPI rank.  Changed the name
of the RAJA ``N_Vector`` library to ``libsundials_nveccudaraja.lib`` from ``libsundials_nvecraja.lib`` to
better reflect that we only support CUDA as a backend for RAJA currently. Several changes were made to
the build system:

-  CMake 3.1.3 is now the minimum required CMake version.

-  Deprecate the behavior of the ``SUNDIALS_INDEX_TYPE`` CMake option and added the ``SUNDIALS_INDEX_SIZE`` CMake option
   to select the ``sunindextype`` integer size.

-  The native CMake FindMPI module is now used to locate an MPI installation.

-  If MPI is enabled and MPI compiler wrappers are not set, the build system will check if ``CMAKE_<language>_COMPILER``
   can compile MPI programs before trying to locate and use an MPI installation.

-  The previous options for setting MPI compiler wrappers and the executable for running MPI programs have been have
   been depreated. The new options that align with those used in native CMake FindMPI module are ``MPI_C_COMPILER``,
   ``MPI_CXX_COMPILER``, ``MPI_Fortran_COMPILER``, and ``MPIEXEC_EXECUTABLE``.

-  When a Fortran name-mangling scheme is needed (e.g., ``ENABLE_LAPACK`` is ``ON``) the build system will infer the
   scheme from the Fortran compiler. If a Fortran compiler is not available or the inferred or default scheme needs to
   be overridden, the advanced options ``SUNDIALS_F77_FUNC_CASE`` and ``SUNDIALS_F77_FUNC_UNDERSCORES`` can be used to
   manually set the name-mangling scheme and bypass trying to infer the scheme.

-  Parts of the main CMakeLists.txt file were moved to new files in the ``src`` and ``example`` directories to make the
   CMake configuration file structure more modular.

Changes in v3.1.2
-----------------

The changes in this minor release include the following:

-  Updated the minimum required version of CMake to 2.8.12 and enabled using rpath by default to locate shared libraries
   on OSX.

-  Fixed Windows specific problem where ``sunindextype`` was not correctly defined when using 64-bit integers for the
   SUNDIALS index type. On Windows ``sunindextype`` is now defined as the MSVC basic type ``__int64``.

-  Added sparse SUNMatrix “Reallocate” routine to allow specification of the nonzero storage.

-  Updated the KLU SUNLinearSolver module to set constants for the two reinitialization types, and fixed a bug in the
   full reinitialization approach where the sparse SUNMatrix pointer would go out of scope on some architectures.

-  Updated the “ScaleAdd” and “ScaleAddI” implementations in the sparse SUNMatrix module to more optimally handle the
   case where the target matrix contained sufficient storage for the sum, but had the wrong sparsity pattern. The sum
   now occurs in-place, by performing the sum backwards in the existing storage. However, it is still more efficient if
   the user-supplied Jacobian routine allocates storage for the sum :math:`I+\gamma J` manually (with zero entries if
   needed).

-  Changed the LICENSE install path to ``instdir/include/sundials``.

Changes in v3.1.1
-----------------

The changes in this minor release include the following:

-  Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers: if “Initialize” was called multiple
   times then the solver memory was reallocated (without being freed).

-  Updated KLU SUNLinearSolver module to use a ``typedef`` for the precision-specific solve function to be used (to
   avoid compiler warnings).

-  Added missing typecasts for some ``(void*)`` pointers (again, to avoid compiler warnings).

-  Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of ``sunindextype`` in one location.

-  Fixed a minor bug in ``KINPrintInfo`` where a case was missing for ``KIN_REPTD_SYSFUNC_ERR`` leading to an undefined
   info message.

-  Added missing ``#include <stdio.h>`` in ``N_Vector`` and ``SUNMatrix`` header files.

-  Fixed an indexing bug in the CUDA ``N_Vector`` implementation of ``N_VWrmsNormMask`` and revised the
   RAJA ``N_Vector`` implementation of ``N_VWrmsNormMask`` to work with mask arrays using values other than zero
   or one. Replaced ``double`` with ``realtype`` in the RAJA vector test functions.

-  Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not require a ``SUNMatrix`` or ``SUNLinearSolver``
   module (e.g., iterative linear solvers or fixed pointer solver).

In addition to the changes above, minor corrections were also made to the example programs, build system, and user
documentation.

Changes in v3.1.0
-----------------

Added ``N_Vector`` print functions that write vector data to a specified file (e.g., ``N_VPrintFile_Serial``).

Added ``make test`` and ``make test_install`` options to the build system for testing SUNDIALS after building with
``make`` and installing with ``make install`` respectively.

Changes in v3.0.0
-----------------

All interfaces to matrix structures and linear solvers have been reworked, and all example programs have been updated.
The goal of the redesign of these interfaces was to provide more encapsulation and ease in the interfacing of custom
linear solvers and interoperability with linear solver libraries. Specific changes include:

-  Added generic SUNMATRIX module with three provided implementations: dense, banded and sparse. These replicate
   previous SUNDIALS Dls and Sls matrix structures in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNMATRIX modules.

-  Added generic ``SUNLinearSolver`` module with eleven provided implementations: SUNDIALS native dense,
   SUNDIALS native banded, LAPACK dense, LAPACK band, KLU, SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, and PCG.
   These replicate previous SUNDIALS generic linear solvers in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNLINEARSOLVER modules.

-  Expanded package-provided direct linear solver (Dls) interfaces and scaled, preconditioned, iterative linear solver
   (Spils) interfaces to utilize generic SUNMATRIX and SUNLINEARSOLVER objects.

-  Removed package-specific, linear solver-specific, solver modules (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since
   their functionality is entirely replicated by the generic Dls/Spils interfaces and SUNLINEARSOLVER/SUNMATRIX modules.
   The exception is CVDIAG, a diagonal approximate Jacobian solver available to CVODE and CVODES.

-  Converted all SUNDIALS example problems to utilize new generic SUNMATRIX and SUNLINEARSOLVER objects, along
   with updated Dls and Spils linear solver interfaces.

-  Added Spils interface routines to ARKode, CVODE, CVODES, IDA and IDAS to allow specification of a user-provided
   "JTSetup" routine. This change supports users who wish to set up data structures for the user-provided
   Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup setup per Newton iteration can be
   amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added – one for CUDA and one for RAJA vectors. These
vectors are supplied to provide very basic support for running on GPU architectures. Users are advised that these
vectors both move all data to the GPU device upon construction, and speedup will only be realized if the user also
conducts the right-hand-side function evaluation on the device. In addition, these vectors assume the problem fits on
one GPU. Further information about RAJA, users are referred to th web site, https://software.llnl.gov/RAJA/. These
additions are accompanied by additions to various interface functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can be configured to be a 32- or 64-bit
integer data index type. ``sunindextype`` is defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``. The Fortran interfaces continue to use ``long int`` for indices,
except for their sparse matrix interface that now uses the new ``sunindextype``. This new flexible capability for index
types includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with either 32-bit or 64-bit capabilities depending how
the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype`` values ``TRUE`` and ``FALSE`` have been
changed to ``SUNTRUE`` and ``SUNFALSE`` respectively.

Temporary vectors were removed from preconditioner setup and solve routines for all packages. It is assumed that all
necessary data for user-provided preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS type information for use in Fortran
programs.

The build system was expanded to support many of the xSDK-compliant keys. The xSDK is a movement in scientific software
to provide a foundation for the rapid and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

Added functions ``SUNDIALSGetVersion`` and ``SUNDIALSGetVersionNumber`` to get SUNDIALS release version
information at runtime.

In addition, numerous changes were made to the build system. These include the addition of separate ``BLAS_ENABLE`` and
``BLAS_LIBRARIES`` CMake variables, additional error checking during CMake configuration, minor bug fixes, and renaming
CMake options to enable/disable examples for greater clarity and an added option to enable/disable Fortran 77 examples.
These changes included changing ``EXAMPLES_ENABLE`` to ``EXAMPLES_ENABLE_C``, changing ``CXX_ENABLE`` to
``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to ``EXAMPLES_ENABLE_F90``, and adding an ``EXAMPLES_ENABLE_F77``
option.

A bug fix was done to correct the fcmix name translation for ``FKIN_SPFGMR``.

Corrections and additions were made to the examples, to installation-related files, and to the user documentation.

Changes in v2.9.0
-----------------

Two additional ``N_Vector`` implementations were added – one for Hypre (parallel) vectors, and one for PETSc vectors.
These additions are accompanied by additions to various interface functions and to user documentation.

Each ``N_Vector`` module now includes a function, ``N_VGetVectorID``, that returns the ``N_Vector`` module name.

The Picard iteration return was chanegd to always return the newest iterate upon success. A minor bug in the line search
was fixed to prevent an infinite loop when the beta condition fails and lamba is below the minimum size.

For each linear solver, the various solver performance counters are now initialized to 0 in both the solver
specification function and in solver ``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

A memory leak was fixed in the banded preconditioner interface. In addition, updates were done to return integers from
linear solver and preconditioner ’free’ functions.

Corrections were made to three Fortran interface functions. The Anderson acceleration scheme was enhanced by use of QR
updating.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant dot product. Various additions and corrections
were made to the interfaces to the sparse solvers KLU and SuperLU_MT, including support for CSR format when using KLU.

The functions FKINCREATE and FKININIT were added to split the FKINMALLOC routine into two pieces. FKINMALLOC remains for
backward compatibility, but documentation for it has been removed.

A new examples was added for use of the OpenMP vector.

Minor corrections and additions were made to the KINSOL solver, to the Fortran interfaces, to the examples, to
installation-related files, and to the user documentation.

Changes in v2.8.0
-----------------

Two major additions were made to the globalization strategy options (``KINSol`` argument ``strategy``). One is
fixed-point iteration, and the other is Picard iteration. Both can be accelerated by use of the Anderson acceleration
method. See the relevant paragraphs in Chapter :numref:`KINSOL.Mathematics`.

Three additions were made to the linear system solvers that are available for use with the KINSOL solver. First,
in the serial case, an interface to the sparse direct solver KLU was added. Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct solver option, to be used with the
serial version of the ``N_Vector`` module. As part of these additions, a sparse matrix (CSC format) structure was added
to KINSOL. Finally, a variation of GMRES called Flexible GMRES was added.

Otherwise, only relatively minor modifications were made to KINSOL:

In function ``KINStop``, two return values were corrected to make the values of ``uu`` and ``fval`` consistent.

A bug involving initialization of ``mxnewtstep`` was fixed. The error affects the case of repeated user calls to
``KINSol`` with no intervening call to ``KINSetMaxNewtonStep``.

A bug in the increments for difference quotient Jacobian approximations was fixed in function ``kinDlsBandDQJac``.

In ``KINLapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to ``smu = mu + ml`` to correct an illegal input
error for ``DGBTRF/DGBTRS``.

In order to avoid possible name conflicts, the mathematical macro and function names ``MIN``, ``MAX``, ``SQR``,
``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and ``RPowerR`` were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``,
``SUNRabs``, ``SUNRsqrt``, ``SUNRexp``, ``SRpowerI``, and ``SUNRpowerR``, respectively. These names occur in both the
solver and in various example programs.

In the FKINSOL module, an incorrect return value ``ier`` in ``FKINfunc`` was fixed.

In the FKINSOL optional input routines ``FKINSETIIN``, ``FKINSETRIN``, and ``FKINSETVIN``, the optional fourth argument
``key_length`` was removed, with hardcoded key string lengths passed to all ``strncmp`` tests.

In all FKINSOL examples, integer declarations were revised so that those which must match a C type ``long int`` are
declared ``INTEGER*8``, and a comment was added about the type match. All other integer declarations are just
``INTEGER``. Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing environments — one for OpenMP, denoted
``NVECTOR_OPENMP``, and one for Pthreads, denoted ``NVECTOR_PTHREADS``.

With this version of SUNDIALS, support and documentation of the Autotools mode of installation is being dropped,
in favor of the CMake mode, which is considered more widely portable.

Changes in v2.7.0
-----------------

One significant design change was made with this release: The problem size and its relatives, bandwidth parameters,
related internal indices, pivot arrays, and the optional output ``lsflag`` have all been changed from type ``int`` to
type ``long int``, except for the problem size and bandwidths in user calls to routines specifying BLAS/LAPACK routines
for the dense/band linear solvers. The function ``NewIntArray`` is replaced by a pair ``NewIntArray``/``NewLintArray``,
for ``int`` and ``long int`` arrays, respectively.

A large number of errors have been fixed. Three major logic bugs were fixed – involving updating the solution vector,
updating the linesearch parameter, and a missing error return. Three minor errors were fixed – involving setting
``etachoice`` in the Matlab/KINSOL interface, a missing error case in ``KINPrintInfo``, and avoiding an
exponential overflow in the evaluation of ``omega``. In each linear solver interface function, the linear solver memory
is freed on an error return, and the ``**Free`` function now includes a line setting to NULL the main memory pointer to
the linear solver memory. In the installation files, we modified the treatment of the macro SUNDIALS_USE_GENERIC_MATH,
so that the parameter GENERIC_MATH_LIB is either defined (with no value) or not defined.

Changes in v2.6.0
-----------------

This release introduces a new linear solver module, based on BLAS and LAPACK for both dense and banded matrices.

The user interface has been further refined. Some of the API changes involve: (a) a reorganization of all linear solver
modules into two families (besides the already present family of scaled preconditioned iterative linear solvers, the
direct solvers, including the new LAPACK-based ones, were also organized into a *direct* family); (b) maintaining a
single pointer to user data, optionally specified through a ``Set``-type function; (c) a general streamlining of the
band-block-diagonal preconditioner module distributed with the solver.

Changes in v2.5.0
-----------------

The main changes in this release involve a rearrangement of the entire SUNDIALS source tree (see
:numref:`KINSOL.Organization`). At the user interface level, the main impact is in the mechanism of including
SUNDIALS header files which must now include the relative path (e.g. ``#include <cvode/cvode.h>``). Additional
changes were made to the build system: all exported header files are now installed in separate subdirectories of the
installation *include* directory.

The functions in the generic dense linear solver (``sundials_dense`` and ``sundials_smalldense``) were modified to work
for rectangular :math:`m \times n` matrices (:math:`m \le n`), while the factorization and solution functions were
renamed to ``DenseGETRF``/``denGETRF`` and ``DenseGETRS``/``denGETRS``, respectively. The factorization and solution
functions in the generic band linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.

Changes in v2.4.0
-----------------

KINSPBCG, KINSPTFQMR, KINDENSE, and KINBAND modules have been added to interface with the Scaled
Preconditioned Bi-CGStab (SPBCG), Scaled Preconditioned Transpose-Free Quasi-Minimal Residual (SPTFQMR),
DENSE, and BAND linear solver modules, respectively. (For details see Chapter :numref:KINSOL.Usage.CC.)
Corresponding additions were made to the Fortran interface module FKINSOL. At the same time, function type names
for Scaled Preconditioned Iterative Linear Solvers were added for the user-supplied Jacobian-times-vector and
preconditioner setup and solve functions.

Regarding the Fortran interface module FKINSOL, optional inputs are now set using ``FKINSETIIN`` (integer inputs),
``FKINSETRIN`` (real inputs), and ``FKINSETVIN`` (vector inputs). Optional outputs are still obtained from the ``IOUT``
and ``ROUT`` arrays which are owned by the user and passed as arguments to ``FKINMALLOC``.

The KINDENSE and KINBAND linear solver modules include support for nonlinear residual monitoring which can
be used to control Jacobian updating.

To reduce the possibility of conflicts, the names of all header files have been changed by adding unique prefixes
(``kinsol_`` and ``sundials_``). When using the default installation procedure, the header files are exported under
various subdirectories of the target ``include`` directory. For more details see Appendix :numref:`Installation`.

Changes in v2.3.0
-----------------

The user interface has been further refined. Several functions used for setting optional inputs were combined into a
single one. Additionally, to resolve potential variable scope issues, all SUNDIALS solvers release user data right
after its use. The build system has been further improved to make it more robust.

Changes in v2.2.1
-----------------

The changes in this minor SUNDIALS release affect only the build system.

Changes in v2.2.0
-----------------

The major changes from the previous version involve a redesign of the user interface across the entire SUNDIALS
suite. We have eliminated the mechanism of providing optional inputs and extracting optional statistics from the solver
through the ``iopt`` and ``ropt`` arrays. Instead, KINSOL now provides a set of routines (with prefix ``KINSet``)
to change the default values for various quantities controlling the solver and a set of extraction routines (with prefix
``KINGet``) to extract statistics after return from the main solver routine. Similarly, each linear solver module
provides its own set of ``Set``- and ``Get``-type routines. For more details see Chapter :numref:KINSOL.Usage.CC.

Additionally, the interfaces to several user-supplied routines (such as those providing Jacobian-vector products and
preconditioner information) were simplified by reducing the number of arguments. The same information that was
previously accessible through such arguments can now be obtained through ``Get``-type functions.

Installation of KINSOL (and all of SUNDIALS) has been completely redesigned and is now based on configure
scripts.

.. _KINSOL.Introduction.reading:

Reading this User Guide
=======================

This user guide is a combination of general usage instructions and specific examples. We expect that some readers will
want to concentrate on the general instructions, while others will refer mostly to the examples, and the organization is
intended to accommodate both styles.

There are different possible levels of usage of KINSOL. The most casual user, with a small nonlinear system, can
get by with reading all of Chapter :numref:`KINSOL.Mathematics`, then Chapter :numref:KINSOL.Usage.CC through
:numref:`KINSOL.Usage.CC` only, and looking at examples in :cite:p:`kinsol_ex`. In a different
direction, a more expert user with a nonlinear system may want to (a) use a package preconditioner
(:numref:`KINSOL.Usage.CC.kin_bbdpre`), (b) supply his/her own Jacobian or preconditioner routines
(:numref:`KINSOL.Usage.CC.user_fct_sim`), (c) supply a new ``N_Vector`` module
(Chapter :numref:`NVectors`), or even (d) supply a different linear solver module
(:numref:`KINSOL.Usage.CC.callable_fct_sim.lin_solv_init` and Chapter :numref:`SUNLinSol`).

The structure of this document is as follows:

-  In Chapter :numref:`KINSOL.Mathematics`, we provide short descriptions of the numerical methods implemented by KINSOL
   for the solution of nonlinear systems.

-  The following chapter describes the structure of the SUNDIALS suite of solvers
   (:numref:`KINSOL.Organization`) and the software organization of the KINSOL solver
   (:numref:`KINSOL.Organization.KINSOL`).

-  Chapter :numref:KINSOL.Usage.CC is the main usage document for KINSOL for C applications. It includes a
   complete description of the user interface for the solution of nonlinear algebraic systems.

-  Chapter :numref:`NVectors` gives a brief overview of the generic ``N_Vector`` module shared among the
   various components of SUNDIALS, and details on the four ``N_Vector`` implementations provided with
   SUNDIALS.

-  Chapter :numref:`SUNMatrix` gives a brief overview of the generic ``SUNMatrix`` module shared among
   the various components of SUNDIALS, and details on the ``SUNMatrix`` implementations provided with
   SUNDIALS.

-  Chapter :numref:`SUNLinSol` gives a brief overview of the generic ``SUNLinearSolver`` module shared among
   the various components of SUNDIALS. This chapter contains details on the ``SUNLinearSolver`` implementations
   provided with SUNDIALS. The chapter also contains details on the ``SUNLinearSolver`` implementations provided with
   SUNDIALS that interface with external linear solver libraries.

-  Finally, in the appendices, we provide detailed instructions for the installation of KINSOL, within the
   structure of SUNDIALS (Appendix :numref:`Installation`), as well as a list of all the constants used for
   input to and output from KINSOL functions (Appendix :numref:`KINSOL.Constants`).

Finally, the reader should be aware of the following notational conventions in this user guide: program listings and
identifiers (such as ``KINInit``) within textual explanations appear in typewriter type style; fields in C
structures (such as *content*) appear in italics; and packages or modules are written in all capitals. Usage and


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.


Acknowledgments
===============

We wish to acknowledge the contributions to previous versions of the KINSOL code and user guide by Allan G.
Taylor.
