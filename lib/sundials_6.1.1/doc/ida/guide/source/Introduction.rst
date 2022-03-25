.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDA.Introduction:

************
Introduction
************

IDA is part of a software family called SUNDIALS: SUite of Nonlinear and
DIfferential/ALgebraic equation Solvers :cite:p:`HBGLSSW:05`.  This suite
consists of CVODE, ARKODE, KINSOL, and IDA, and variants of these with
sensitivity analysis capabilities, CVODES and IDAS.

IDA is a general purpose solver for the initial value problem (IVP) for systems
of differential-algebraic equations (DAEs). The name IDA stands for Implicit
Differential-Algebraic solver. IDA is based on DASPK :cite:p:`BHP:94,BHP:98`,
but is written in ANSI-standard C rather than Fortran77.  Its most notable
features are that,
(1) in the solution of the underlying nonlinear system at each time step, it
offers a choice of Newton/direct methods and a choice of Inexact Newton/Krylov
(iterative) methods; and
(2) it is written in a *data-independent* manner in that it acts on generic
vectors and matrices without any assumptions on the underlying organization of
the data.  Thus IDA shares significant modules previously written within CASC at
LLNL to support the ordinary differential equation (ODE) solvers CVODE
:cite:p:`cvode_ug,CoHi:96` and PVODE :cite:p:`ByHi:98,ByHi:99`, and also the
nonlinear system solver KINSOL :cite:p:`kinsol_ug`.

At present, IDA may utilize a variety of Krylov methods provided in SUNDIALS
that can be used in conjuction with Newton iteration: these include the GMRES
(Generalized Minimal RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized
Minimum RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate Gradient Stabilized)
:cite:p:`Van:92`, TFQMR (Transpose-Free Quasi-Minimal Residual)
:cite:p:`Fre:93`, and PCG (Preconditioned Conjugate Gradient) :cite:p:`HeSt:52`
linear iterative methods. As Krylov methods, these require little matrix storage
for solving the Newton equations as compared to direct methods. However, the
algorithms allow for a user-supplied preconditioner, and, for most
problems, preconditioning is essential for an efficient solution.

For very large DAE systems, the Krylov methods are preferable over direct linear
solver methods, and are often the only feasible choice.  Among the Krylov
methods in SUNDIALS, we recommend GMRES as the best overall choice. However,
users are encouraged to compare all options, especially if encountering
convergence failures with GMRES.  Bi-CGFStab and TFQMR have an advantage in
storage requirements, in that the number of workspace vectors they require is
fixed, while that number for GMRES depends on the desired Krylov subspace
size. FGMRES has an advantage in that it is designed to support preconditioners
that vary between iterations (e.g. iterative methods). PCG exhibits rapid
convergence and minimal workspace vectors, but only works for symmetric linear
systems.

..
   There are several motivations for choosing the C language for IDA.  First, a
   general movement away from Fortran and toward C in scientific computing was
   apparent. Second, the pointer, structure, and dynamic memory allocation features
   in C are extremely useful in software of this complexity, with the great variety
   of method options offered.  Finally, we prefer C over C++ for IDA because of the
   wider availability of C compilers, the potentially greater efficiency of C, and
   the greater ease of interfacing the solver to applications written in extended
   Fortran.

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

**IDA**

The IDA Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
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
| ``IDASpilsSetLinearSolver``     | ``IDASetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``IDASpilsSetPreconditioner``   | ``IDASetPreconditioner``       |
+---------------------------------+--------------------------------+
| ``IDASpilsSetJacTimes``         | ``IDASetJacTimes``             |
+---------------------------------+--------------------------------+
| ``IDASpilsSetEpsLin``           | ``IDASetEpsLin``               |
+---------------------------------+--------------------------------+
| ``IDASpilsSetIncrementFactor``  | ``IDASetIncrementFactor``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetWorkSpace``        | ``IDAGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumPrecEvals``     | ``IDAGetNumPrecEvals``         |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumPrecSolves``    | ``IDAGetNumPrecSolves``        |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumLinIters``      | ``IDAGetNumLinIters``          |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumConvFails``     | ``IDAGetNumLinConvFails``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumJTSetupEvals``  | ``IDAGetNumJTSetupEvals``      |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumJtimesEvals``   | ``IDAGetNumJtimesEvals``       |
+---------------------------------+--------------------------------+
| ``IDASpilsGetNumResEvals``      | ``IDAGetNumLinResEvals``       |
+---------------------------------+--------------------------------+
| ``IDASpilsGetLastFlag``         | ``IDAGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``IDASpilsGetReturnFlagName``   | ``IDAGetLinReturnFlagName``    |
+---------------------------------+--------------------------------+
| ``IDADlsSetLinearSolver``       | ``IDASetLinearSolver``         |
+---------------------------------+--------------------------------+
| ``IDADlsSetJacFn``              | ``IDASetJacFn``                |
+---------------------------------+--------------------------------+
| ``IDADlsGetWorkSpace``          | ``IDAGetLinWorkSpace``         |
+---------------------------------+--------------------------------+
| ``IDADlsGetNumJacEvals``        | ``IDAGetNumJacEvals``          |
+---------------------------------+--------------------------------+
| ``IDADlsGetNumResEvals``        | ``IDAGetNumLinResEvals``       |
+---------------------------------+--------------------------------+
| ``IDADlsGetLastFlag``           | ``IDAGetLastLinFlag``          |
+---------------------------------+--------------------------------+
| ``IDADlsGetReturnFlagName``     | ``IDAGetLinReturnFlagName``    |
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

The :ref:`RAJA N_Vector <NVectors.RAJA>` implementation has been updated to
support the SYCL backend in addition to the CUDA and HIP backends. Users can
choose the backend when configuring SUNDIALS by using the
:cmakeop:`SUNDIALS_RAJA_BACKENDS` CMake variable. This module remains
experimental and is subject to change from version to version.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix
and the linear solver support general dense linear systems as well as block
diagonal linear systems. See :numref:`SUNLinSol.OneMklDense` for more
details. This module is experimental and is subject to change from version to
version.

Added a new *optional* function to the ``SUNLinearSolver`` API,
:c:func:`SUNLinSolSetZeroGuess`, to indicate that the next call to
:c:func:`SUNLinSolSolve` will be made with a zero initial guess.
``SUNLinearSolver`` implementations that do not use the
:c:func:`SUNLinSolNewEmpty` constructor will, at a minimum, need set the
``setzeroguess`` function pointer in the linear solver ``ops`` structure to
``NULL``. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

IDA now supports a new "matrix-embedded" ``SUNLinearSolver`` type. This type
supports user-supplied ``SUNLinearSolver`` implementations that set up and solve
the specified linear system at each linear solve call. Any matrix-related data
structures are held internally to the linear solver itself, and are not provided
by the SUNDIALS package.

Added the function :c:func:`IDASetNlsResFn` to supply an alternative residual
side function for use within nonlinear system function evaluations.

The installed ``SUNDIALSConfig.cmake`` file now supports the ``COMPONENTS``
option to ``find_package``.

A bug was fixed in :c:func:`SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

Changes in v5.7.0
-----------------

A new ``N_Vector`` implementation based on the SYCL abstraction layer has been
added targeting Intel GPUs. At present the only SYCL compiler supported is the
DPC++ (Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details.
This module is considered experimental and is subject to major changes even in
minor releases.

A new ``SUNMatrix`` and ``SUNLinearSolver`` implementation were added to
interface with the MAGMA linear algebra library. Both the matrix and the linear
solver support general dense linear systems as well as block diagonal linear
systems, and both are targeted at GPUs (AMD or NVIDIA). See
:numref:`SUNLinSol.MagmaDense` for more details.

Changes in v5.6.1
-----------------

Fixed a bug in the SUNDIALS CMake which caused an error if the
:cmakeop:`CMAKE_CXX_STANDARD` and :cmakeop:`SUNDIALS_RAJA_BACKENDS` options were
not provided.

Fixed some compiler warnings when using the IBM XL compilers.

Changes in v5.6.0
-----------------

A new ``N_Vector`` implementation based on the AMD ROCm HIP platform has been
added. This vector can target NVIDIA or AMD GPUs. See :numref:`NVectors.Hip` for
more details. This module is considered experimental and is subject to change
from version to version.

The :ref:`NVECTOR_RAJA <NVectors.RAJA>` implementation has been updated to
support the HIP backend in addition to the CUDA backend. Users can choose the
backend when configuring SUNDIALS by using the :cmakeop:`SUNDIALS_RAJA_BACKENDS`
CMake variable. This module remains experimental and is subject to change from
version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the
``N_Vector`` API. This operation is useful for :c:type:`N_Vectors` that utilize
dual memory spaces, e.g. the native SUNDIALS CUDA ``N_Vector``.

The :ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` and
:ref:`SUNLINEARSOLVER_CUSOLVERSP_BATCHQR <SUNLinSol.cuSolverSp>` implementations
no longer require the SUNDIALS CUDA ``N_Vector``. Instead, they require that the
vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and
that the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA
device pointer.

Changes in v5.5.0
-----------------

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes should be
fully backwards compatible for almost all users. SUNDIALS now exports CMake
targets and installs a ``SUNDIALSConfig.cmake`` file.

Added support for SuperLU_DIST 6.3.0 or newer.

Changes in v5.4.0
-----------------

Added the function :c:func:`IDASetLSNormFactor` to specify the factor for
converting between integrator tolerances (WRMS norm) and linear solver
tolerances (L2 norm) i.e., ``tol_L2 = nrmfac * tol_WRMS``.

The expected behavior of :c:func:`SUNNonlinSolGetNumIters` and
:c:func:`SUNNonlinSolGetNumConvFails` in the ``SUNNonlinearSolver`` API have
been updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided ``SUNNonlinearSolver``
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retreived by calling
:c:func:`IDAGetNumNonlinSolvIters`, the cumulative number of failures with
:c:func:`IDAGetNumNonlinSolvConvFails`, or both with
:c:func:`IDAGetNonlinSolvStats`.

A new API, ``SUNMemoryHelper``, was added to support **GPU users** who have
complex memory management needs such as using memory pools. This is paired with
new constructors for the :ref:`NVECTOR_CUDA <NVectors.CUDA>` and
:ref:`NVECTOR_RAJA <NVectors.RAJA>` modules that accept a ``SUNMemoryHelper``
object. Refer to :numref:`SUNDIALS.GPU` and :numref:`SUNMemory` for more
information.

The :ref:`NVECTOR_RAJA <NVectors.RAJA>` module has been updated to mirror the
:ref:`NVECTOR_CUDA <NVectors.CUDA>` module.  Notably, the update adds managed
memory support to the :ref:`NVECTOR_RAJA <NVectors.RAJA>` module.  Users of the
module will need to update any calls to the :c:func:`N_VMake_Raja` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

The :ref:`NVECTOR_TRILINOS <NVectors.NVTrilinos>` module has been updated to
work with Trilinos 12.18+. This update changes the local ordinal type to always
be an ``int``.

Added support for CUDA v11.

Changes in v5.3.0
-----------------

Fixed a bug in the iterative linear solver modules where an error is not
returned if the ATimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added a new function :c:func:`IDAGetNonlinearSystemData` which advanced users
might find useful if providing a custom :c:type:`SUNNonlinSolSysFn`.

Added the ability to control the CUDA kernel launch parameters for the
:ref:`NVECTOR_CUDA <NVectors.CUDA>` and
:ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` modules. These modules remain
experimental and are subject to change from version to version.  In addition,
the :ref:`NVECTOR_CUDA <NVectors.CUDA>` kernels were rewritten to be more
flexible. Most users should see equivalent performance or some improvement, but
a select few may observe minor performance degradation with the default
settings. Users are encouraged to contact the SUNDIALS team about any
performance changes that they notice.

Added new capabilities for monitoring the solve phase in the
:ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
and :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>` modules, and the
SUNDIALS iterative linear solver modules. SUNDIALS must be built with the CMake
option :cmakeop:`SUNDIALS_BUILD_WITH_MONITORING` to use these capabilities.

Added the optional function :c:func:`IDASetJacTimesResFn` to specify an
alternative residual function for computing Jacobian-vector products with the
internal difference quotient approximation.

Changes in v5.2.0
-----------------

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set :cmakeop:`CMAKE_Fortran_COMPILER` to ``f2003``,
``xlf2003``, or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attributes missing on some SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, :ref:`SUNMATRIX_CUSPARSE
<SUNMatrix.cuSparse>`, that interfaces to the sparse matrix implementation from
the NVIDIA cuSPARSE library. In addition, the :ref:`SUNLINSOL_CUSOLVER_BATCHQR
<SUNLinSol.cuSolverSp>` linear solver has been updated to use this matrix,
therefore, users of this module will need to update their code.  These modules
are still considered to be experimental, thus they are subject to breaking
changes even in minor releases.

The function :c:func:`IDASetLinearSolutionScaling` was added to enable or
disable the scaling applied to linear system solutions with matrix-based linear
solvers to account for a lagged value of :math:`\alpha` in the linear system
matrix :math:`J = \frac{\partial F}{\partial y} + \alpha\frac{\partial
F}{\partial \dot{y}}`.  Scaling is enabled by default when using a matrix-based
linear solver.

Changes in v5.1.0
-----------------

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake variables
:cmakeop:`PETSC_INCLUDES` and :cmakeop:`PETSC_LIBRARIES` instead of
:cmakeop:`PETSC_DIR`.

Added a new build system option, :cmakeop:`CUDA_ARCH`, that can be used to
specify the CUDA architecture to compile for.

Added two utility functions, :f:func:`FSUNDIALSFileOpen` and
:f:subr:`FSUNDIALSFileClose` for creating/destroying file pointers that are
useful when using the Fortran 2003 interfaces.

Changes in v5.0.0
-----------------

Build system changes
^^^^^^^^^^^^^^^^^^^^

* Increased the minimum required CMake version to 3.5 for most SUNDIALS
  configurations, and 3.10 when CUDA or OpenMP with device offloading are
  enabled.

* The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been
  removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
  third party libraries that require linking to BLAS, the path to the BLAS
  library should be included in the ``*_LIBRARIES`` variable for the third party
  library *e.g.*, :cmakeop:`SUPERLUDIST_LIBRARIES` when enabling SuperLU_DIST.

* Fixed a bug in the build system that prevented the
  :ref:`NVECTOR_PTHREADS <NVectors.Pthreads>` module from being built.

NVECTOR module changes
^^^^^^^^^^^^^^^^^^^^^^

* Two new functions were added to aid in creating custom ``N_Vector``
  objects. The constructor :c:func:`N_VNewEmpty` allocates an "empty" generic
  ``N_Vector`` with the object’s content pointer and the function pointers in
  the operations structure initialized to ``NULL``. When used in the constructor
  for custom objects this function will ease the introduction of any new
  optional operations to the ``N_Vector`` API by ensuring only required
  operations need to be set.  Additionally, the function :c:func:`N_VCopyOps`
  has been added to copy the operation function pointers between vector
  objects. When used in clone routines for custom vector objects these functions
  also will ease the introduction of any new optional operations to the
  ``N_Vector`` API by ensuring all operations are copied when cloning
  objects. See :numref:`NVectors.Description.utilities` for more details.

* Two new ``N_Vector`` implementations,
  :ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` and
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>`, have been created to
  support flexible partitioning of solution data among different processing
  elements (e.g., CPU + GPU) or for multi-physics problems that couple distinct
  MPI-based simulations together. This implementation is accompanied by
  additions to user documentation and SUNDIALS examples. See
  :numref:`NVectors.ManyVector` and :numref:`NVectors.MPIManyVector` for more
  details.

* One new required vector operation and ten new optional vector operations have
  been added to the ``N_Vector`` API. The new required operation,
  :c:func:`N_VGetLength`, returns the global length of an ``N_Vector``.  The
  optional operations have been added to support the new
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` implementation. The
  operation :c:func:`N_VGetCommunicator` must be implemented by subvectors that
  are combined to create an
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>`, but is not used outside
  of this context. The remaining nine operations are optional local reduction
  operations intended to eliminate unnecessary latency when performing vector
  reduction operations (norms, etc.) on distributed memory systems. The optional
  local reduction vector operations are :c:func:`N_VDotProdLocal`,
  :c:func:`N_VMaxNormLocal`, :c:func:`N_VMinLocal`, :c:func:`N_VL1NormLocal`,
  :c:func:`N_VWSqrSumLocal`, :c:func:`N_VWSqrSumMaskLocal`,
  :c:func:`N_VInvTestLocal`, :c:func:`N_VConstrMaskLocal`, and
  :c:func:`N_VMinQuotientLocal`. If an ``N_Vector`` implementation defines any
  of the local operations as ``NULL``, then the
  :ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` will call standard
  ``N_Vector`` operations to complete the computation. See
  :numref:`NVectors.Ops.Local` for more details.

* An additional ``N_Vector`` implementation, :ref:`NVECTOR_MPIPLUSX
  <NVectors.MPIPlusX>`, has been created to support the MPI+X paradigm where X
  is a type of on-node parallelism (*e.g.*, OpenMP, CUDA). The implementation is
  accompanied by additions to user documentation and SUNDIALS examples. See
  :numref:`NVectors.MPIPlusX` for more details.

* The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` and :ref:`NVECTOR_RAJA <NVectors.RAJA>`
  implementations respectively. Accordingly, the ``nvector_mpicuda.h``,
  ``nvector_mpiraja.h``, ``libsundials_nvecmpicuda.lib``, and
  ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should use
  the :ref:`NVECTOR_MPIPLUSX <NVectors.MPIPlusX>` module coupled in conjunction
  with the :ref:`NVECTOR_CUDA <NVectors.CUDA>` or :ref:`NVECTOR_RAJA
  <NVectors.RAJA>` modules to replace the functionality. The necessary changes
  are minimal and should require few code modifications. See the programs in
  ``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to
  use the :ref:`NVECTOR_MPIPLUSX <NVectors.MPIPlusX>` module with the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` and :ref:`NVECTOR_RAJA <NVectors.RAJA>`
  modules respectively.

* Fixed a memory leak in the :ref:`NVECTOR_PETSC <NVectors.NVPETSc>` module
  clone function.

* Made performance improvements to the :ref:`NVECTOR_CUDA <NVectors.CUDA>`
  module. Users who utilize a non-default stream should no longer see default
  stream synchronizations after memory transfers.

* Added a new constructor to the :ref:`NVECTOR_CUDA <NVectors.CUDA>` module that
  allows a user to provide custom allocate and free functions for the vector
  data array and internal reduction buffer. See :numref:`NVectors.CUDA` for more
  details.

* Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See
  :numref:`NVectors` for more details on how to use the interfaces.

* Added three new ``N_Vector`` utility functions,
  :c:func:`FN_VGetVecAtIndexVectorArray`,
  :c:func:`FN_VSetVecAtIndexVectorArray`, and :c:func:`FN_VNewVectorArray`, for
  working with ``N_Vector`` arrays when using the Fortran 2003 interfaces.  See
  :numref:`NVectors.Description.utilities` for more details.

SUNMatrix module changes
^^^^^^^^^^^^^^^^^^^^^^^^

* Two new functions were added to aid in creating custom ``SUNMatrix``
  objects. The constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic
  ``SUNMatrix`` with the object’s content pointer and the function pointers in
  the operations structure initialized to ``NULL``. When used in the constructor
  for custom objects this function will ease the introduction of any new
  optional operations to the ``SUNMatrix`` API by ensuring only required
  operations need to be set.  Additionally, the function :c:func:`SUNMatCopyOps`
  has been added to copy the operation function pointers between matrix
  objects. When used in clone routines for custom matrix objects these functions
  also will ease the introduction of any new optional operations to the
  ``SUNMatrix`` API by ensuring all operations are copied when cloning
  objects. See :numref:`SUNMatrix.Description` for more details.

* A new operation, :c:func:`SUNMatMatvecSetup`, was added to the ``SUNMatrix``
  API to perform any setup necessary for computing a matrix-vector product. This
  operation is useful for ``SUNMatrix`` implementations which need to prepare
  the matrix itself, or communication structures before performing the
  matrix-vector product. Users who have implemented custom ``SUNMatrix`` modules
  will need to at least update their code to set the corresponding ``ops``
  structure member, ``matvecsetup``, to ``NULL``. See
  :numref:`SUNMatrix.Description` for more details.

* The generic ``SUNMatrix`` API now defines error codes to be returned by
  ``SUNMatrix`` operations. Operations which return an integer flag indicating
  success/failure may return different values than previously.

* A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to
  facilitate the use of the SuperLU_DIST library with SUNDIALS. See
  :numref:`SUNMatrix.SLUNRloc` for more details.

* Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See
  :numref:`SUNMatrix` for more details on how to use the interfaces.

SUNLinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* A new function was added to aid in creating custom ``SUNLinearSolver``
  objects.  The constructor :c:func:`SUNLinSolNewEmpty` allocates an "empty"
  generic ``SUNLinearSolver`` with the object’s content pointer and the function
  pointers in the operations structure initialized to ``NULL``. When used in the
  constructor for custom objects this function will ease the introduction of any
  new optional operations to the ``SUNLinearSolver`` API by ensuring only
  required operations need to be set. See :numref:`SUNLinSol.API.Custom` for
  more details.

* The return type of the ``SUNLinearSolver`` API function
  :c:func:`SUNLinSolLastFlag` has changed from ``long int`` to ``sunindextype``
  to be consistent with the type used to store row indices in dense and banded
  linear solver modules.

* Added a new optional operation to the ``SUNLinearSolver`` API,
  :c:func:`SUNLinSolGetID`, that returns a ``SUNLinearSolver_ID`` for
  identifying the linear solver module.

* The ``SUNLinearSolver`` API has been updated to make the initialize and setup
  functions optional.

* A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to
  facilitate the use of the SuperLU_DIST library with SUNDIALS. See
  :numref:`SUNLinSol.SuperLUDIST` for more details.

* Added a new ``SUNLinearSolver`` implementation,
  SUNLinearSolver_cuSolverSp_batchQR, which leverages the NVIDIA cuSOLVER sparse
  batched QR method for efficiently solving block diagonal linear systems on
  NVIDIA GPUs. See :numref:`SUNLinSol.cuSolverSp` for more details.

* Added three new accessor functions to the SUNLINSOL_KLU module,
  :c:func:`SUNLinSol_KLUGetSymbolic`, :c:func:`SUNLinSol_KLUGetNumeric`, and
  :c:func:`SUNLinSol_KLUGetCommon`, to provide user access to the underlying KLU
  solver structures. See :numref:`SUNLinSol.KLU` for more details.

* Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules.  See
  :numref:`SUNLinSol` for more details on how to use the interfaces.

SUNNonlinearSolver module changes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* A new function was added to aid in creating custom ``SUNNonlinearSolver``
  objects. The constructor :c:func:`SUNNonlinSolNewEmpty` allocates an "empty"
  generic ``SUNNonlinearSolver`` with the object’s content pointer and the
  function pointers in the operations structure initialized to ``NULL``. When
  used in the constructor for custom objects this function will ease the
  introduction of any new optional operations to the ``SUNNonlinearSolver`` API
  by ensuring only required operations need to be set. See
  :numref:`SUNNonlinSol.API.Custom` for more details.

* To facilitate the use of user supplied nonlinear solver convergence test
  functions the :c:type:`SUNNonlinSolSetConvTestFn` function in the
  ``SUNNonlinearSolver`` API has been updated to take a ``void*`` data pointer
  as input. The supplied data pointer will be passed to the nonlinear solver
  convergence test function on each call.

* The inputs values passed to the first two inputs of the
  :c:func:`SUNNonlinSolSolve` function in the ``SUNNonlinearSolver`` have been
  changed to be the predicted state and the initial guess for the correction to
  that state. Additionally, the definitions of :c:type:`SUNNonlinSolLSetupFn`
  and :c:type:`SUNNonlinSolLSolveFn` in the ``SUNNonlinearSolver`` API have been
  updated to remove unused input parameters. For more information see
  :numref:`SUNNonlinSol`.

* Added a new ``SUNNonlinearSolver`` implementation,
  :ref:`SUNNONLINSOL_PETSC <SUNNonlinSol.PetscSNES>`, which interfaces to the
  PETSc SNES nonlinear solver API. See :numref:`SUNNonlinSol.PetscSNES` for more
  details.

* Added new Fortran 2003 interfaces for most ``SUNNonlinearSolver`` modules. See
  :numref:`SUNNonlinSol` for more details on how to use the interfaces.

IDA changes
^^^^^^^^^^^

* A bug was fixed in the IDA linear solver interface where an incorrect
  Jacobian-vector product increment was used with iterative solvers other than
  :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>` and
  :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>`.

* Fixed a memeory leak in FIDA when not using the default nonlinear solver.

* Removed extraneous calls to :c:func:`N_VMin` for simulations where the scalar
  valued absolute tolerance, or all entries of the vector-valued absolute
  tolerance array, are strictly positive. In this scenario, IDA will remove at
  least one global reduction per time step.

* The IDALS interface has been updated to only zero the Jacobian matrix before
  calling a user-supplied Jacobian evaluation function when the attached linear
  solver has type ``SUNLINEARSOLVER_DIRECT``.

* Added the new functions, :c:func:`IDAGetCurrentCj`, :c:func:`IDAGetCurrentY`,
  :c:func:`IDAGetCurrentYp`, :c:func:`IDAComputeY`, and :c:func:`IDAComputeYp`
  which may be useful to users who choose to provide their own nonlinear solver
  implementations.

* Added a Fortran 2003 interface to IDA. See :numref:`SUNDIALS.Fortran` for more
  details.

Changes in v4.1.0
-----------------

An additional ``N_Vector`` implementation was added for the TPETRA vector from
the TRILINOS library to facilitate interoperability between SUNDIALS and
TRILINOS. This implementation is accompanied by additions to user documentation
and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option
:cmakeop:`EXAMPLES_ENABLE_CUDA` enables all examples that use CUDA including the
RAJA examples with a CUDA back end (if the RAJA ``N_Vector`` is enabled).

The implementation header file ``ida_impl.h`` is no longer installed. This means
users who are directly manipulating the ``IDAMem`` structure will need to update
their code to use IDA’s public API.

Python is no longer required to run ``make test`` and ``make test_install``.

Changes in v4.0.2
-----------------

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source
file. The symbols are now included in the IDA library, ``libsundials_ida``.

Changes in v4.0.1
-----------------

No changes were made in this release.

Changes in v4.0.0
-----------------

IDA’s previous direct and iterative linear solver interfaces, IDADLS and
IDASPILS, have been merged into a single unified linear solver interface, IDALS,
to support any valid ``SUNLinearSolver`` module.  This includes the "DIRECT" and
"ITERATIVE" types as well as the new "MATRIX_ITERATIVE" type. Details regarding
how IDALS utilizes linear solvers of each type as well as discussion regarding
intended use cases for user-supplied ``SUNLinearSolver`` implementations are
included in :numref:`SUNLinSol`. All IDA example programs and the standalone
linear solver examples have been updated to use the unified linear solver
interface.

The unified interface for the new IDALS module is very similar to the previous
IDADLS and IDASPILS interfaces. To minimize challenges in user migration to the
new names, the previous C and Fortran routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. Additionally, we note that Fortran users, however, may need to
enlarge their ``iout`` array of optional integer outputs, and update the indices
that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver``
implementations have been updated to follow the naming convention ``SUNLinSol_``
where ``*`` is the name of the linear solver. The new names are
:c:func:`SUNLinSol_Band`, :c:func:`SUNLinSol_Dense`, :c:func:`SUNLinSol_KLU`,
:c:func:`SUNLinSol_LapackBand`, :c:func:`SUNLinSol_LapackDense`,
:c:func:`SUNLinSol_PCG`, :c:func:`SUNLinSol_SPBCGS`, :c:func:`SUNLinSol_SPFGMR`,
:c:func:`SUNLinSol_SPGMR`, :c:func:`SUNLinSol_SPTFQMR`, and
:c:func:`SUNLinSol_SuperLUMT`. Solver-specific "set" routine names have been
similarly standardized. To minimize challenges in user migration to the new
names, the previous routine names may still be used; these will be deprecated in
future releases, so we recommend that users migrate to the new names soon. All
IDA example programs and the standalone linear solver examples have been updated
to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the storage
upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the ``SUNNonlinearSolver`` API. This API will ease the
addition of new nonlinear solver options and allow for external or user-supplied
nonlinear solvers. The ``SUNNonlinearSolver`` API and SUNDIALS provided modules
are described in :numref:`SUNNonlinSol` and follow the same object oriented
design and implementation used by the ``N_Vector``, ``SUNMatrix``, and
``SUNLinearSolver`` modules. Currently two ``SUNNonlinearSolver``
implementations are provided, :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
and
:ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`. These replicate the
previous integrator specific implementations of a Newton iteration and a
fixed-point iteration (previously referred to as a functional iteration),
respectively. Note the :ref:`SUNNONLINSOL_FIXEDPOINT <SUNNonlinSol.FixedPoint>`
module can optionally utilize Anderson’s method to accelerate
convergence. Example programs using each of these nonlinear solver modules in a
standalone manner have been added and all IDA example programs have been updated
to use generic ``SUNNonlinearSolver`` modules.

By default IDA uses the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>`
module. Since IDA previously only used an internal implementation of a Newton
iteration no changes are required to user programs and functions for setting the
nonlinear solver options (e.g., :c:func:`IDASetMaxNonlinIters`) or getting
nonlinear solver statistics (e.g., :c:func:`IDAGetNumNonlinSolvIters`) remain
unchanged and internally call generic ``SUNNonlinearSolver`` functions as
needed. While SUNDIALS includes a fixed-point nonlinear solver module, it is not
currently supported in IDA. For details on attaching a user-supplied nonlinear
solver to IDA see :numref:IDA.Usage.CC. Additionally, the example program
``idaRoberts_dns.c`` explicitly creates an attaches a :ref:`SUNNONLINSOL_NEWTON
<SUNNonlinSol.Newton>` object to demonstrate the process of creating and
attaching a nonlinear solver module (note this is not necessary in general as
IDA uses the :ref:`SUNNONLINSOL_NEWTON <SUNNonlinSol.Newton>` module by
default).

Three fused vector operations and seven vector array operations have been added
to the ``N_Vector`` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating an
``N_Vector`` (see :numref:`NVectors` for more details). The new operations are
intended to increase data reuse in vector operations, reduce parallel
communication on distributed memory systems, and lower the number of kernel
launches on systems with accelerators. The fused operations are
:c:func:`N_VLinearCombination`, :c:func:`N_VScaleAddMulti`, and
:c:func:`N_VDotProdMulti` and the vector array operations are
:c:func:`N_VLinearCombinationVectorArray`, :c:func:`N_VScaleVectorArray`,
:c:func:`N_VConstVectorArray`, :c:func:`N_VWrmsNormVectorArray`,
:c:func:`N_VWrmsNormMaskVectorArray`, :c:func:`N_VScaleAddMultiVectorArray`, and
:c:func:`N_VLinearCombinationVectorArray`.

If an ``N_Vector`` implementation defines any of these operations as ``NULL``,
then standard ``N_Vector`` operations will automatically be called as necessary
to complete the computation.

Multiple updates to :ref:`NVECTOR_CUDA <NVectors.CUDA>` were made:

* Changed :c:func:`N_VGetLength_Cuda` to return the global vector length instead
  of the local vector length.

* Added :c:func:`N_VGetLocalLength_Cuda` to return the local vector length.

* Added :c:func:`N_VGetMPIComm_Cuda` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

* Changed the :c:func:`N_VMake_Cuda` function to take a host data pointer and a
  device data pointer instead of an ``N_VectorContent_Cuda`` object.

* Added the ability to set the ``cudaStream_t`` used for execution of the
  :ref:`NVECTOR_CUDA <NVectors.CUDA>` kernels. See the function
  :c:func:`N_VSetCudaStreams_Cuda`.

* Added :c:func:`N_VNewManaged_Cuda`, :c:func:`N_VMakeManaged_Cuda`, and
  :c:func:`N_VIsManagedMemory_Cuda` functions to accommodate using managed
  memory with the :ref:`NVECTOR_CUDA <NVectors.CUDA>`.

Multiple changes to :ref:`NVECTOR_RAJA <NVectors.RAJA>` were made:

* Changed :c:func:`N_VGetLength_Raja` to return the global vector length instead
  of the local vector length.

* Added :c:func:`N_VGetLocalLength_Raja` to return the local vector length.

* Added :c:func:`N_VGetMPIComm_Raja` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device offloading
has been added, :ref:`NVECTOR_OPENMPDEV <NVectors.OpenMPDEV>`. See
:numref:`NVectors.OpenMPDEV` for more details.

Changes in v3.2.1
-----------------

The changes in this minor release include the following:

* Fixed a bug in the :ref:`CUDA N_Vector <NVectors.CUDA>` where the
  :c:func:`N_VInvTest` operation could write beyond the allocated vector data.

* Fixed library installation path for multiarch systems. This fix changes the
  default library installation path to
  ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from
  ``CMAKE_INSTALL_PREFIX/lib``. Note :cmakeop:`CMAKE_INSTALL_LIBDIR` is
  automatically set, but is available as a CMake option that can be modified.

Changes in v3.2.0
-----------------

Fixed a problem with setting ``sunindextype`` which would occur with some
compilers (e.g. armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA ``N_Vector`` library to
``libsundials_nveccudaraja.lib`` from ``libsundials_nvecraja.lib`` to better
reflect that we only support CUDA as a backend for RAJA currently.

Several changes were made to the build system:

* CMake 3.1.3 is now the minimum required CMake version.

* Deprecate the behavior of the :cmakeop:`SUNDIALS_INDEX_TYPE` CMake option and
  added the :cmakeop:`SUNDIALS_INDEX_SIZE` CMake option to select the
  ``sunindextype`` integer size.

* The native CMake FindMPI module is now used to locate an MPI installation.

* If MPI is enabled and MPI compiler wrappers are not set, the build system will
  check if ``CMAKE_<language>_COMPILER`` can compile MPI programs before trying
  to locate and use an MPI installation.

* The previous options for setting MPI compiler wrappers and the executable for
  running MPI programs have been have been depreated. The new options that align
  with those used in native CMake FindMPI module are :cmakeop:`MPI_C_COMPILER`,
  :cmakeop:`MPI_CXX_COMPILER`, :cmakeop:`MPI_Fortran_COMPILER`, and
  :cmakeop:`MPIEXEC_EXECUTABLE`.

* When a Fortran name-mangling scheme is needed (e.g., :cmakeop:`ENABLE_LAPACK`
  is ``ON``) the build system will infer the scheme from the Fortran compiler.
  If a Fortran compiler is not available or the inferred or default scheme needs
  to be overridden, the advanced options :cmakeop:`SUNDIALS_F77_FUNC_CASE` and
  :cmakeop:`SUNDIALS_F77_FUNC_UNDERSCORES` can be used to manually set the
  name-mangling scheme and bypass trying to infer the scheme.

* Parts of the main CMakeLists.txt file were moved to new files in the ``src``
  and ``example`` directories to make the CMake configuration file structure
  more modular.

Changes in v3.1.2
-----------------

The changes in this minor release include the following:

* Updated the minimum required version of CMake to 2.8.12 and enabled using
  rpath by default to locate shared libraries on OSX.

* Fixed Windows specific problem where ``sunindextype`` was not correctly
  defined when using 64-bit integers for the SUNDIALS index type. On Windows
  ``sunindextype`` is now defined as the MSVC basic type ``__int64``.

* Added sparse SUNMatrix "Reallocate" routine to allow specification of the
  nonzero storage.

* Updated the KLU SUNLinearSolver module to set constants for the two
  reinitialization types, and fixed a bug in the full reinitialization approach
  where the sparse SUNMatrix pointer would go out of scope on some
  architectures.

* Updated the :c:func:`SUNMatScaleAdd` and :c:func:`SUNMatScaleAddI`
  implementations in the sparse SUNMatrix module to more optimally handle the
  case where the target matrix contained sufficient storage for the sum, but had
  the wrong sparsity pattern. The sum now occurs in-place, by performing the sum
  backwards in the existing storage. However, it is still more efficient if the
  user-supplied Jacobian routine allocates storage for the sum
  :math:`I+\gamma J` manually (with zero entries if needed).

* Changed the LICENSE install path to ``instdir/include/sundials``.

Changes in v3.1.1
-----------------

The changes in this minor release include the following:

* Fixed a potential memory leak in the :ref:`SUNLINSOL_SPGMR <SUNLinSol.SPGMR>`
  and :ref:`SUNLINSOL_SPFGMR <SUNLinSol.SPFGMR>` linear solvers: if
  "Initialize" was called multiple times then the solver memory was reallocated
  (without being freed).

* Updated KLU ``SUNLinearSolver`` module to use a ``typedef`` for the
  precision-specific solve function to be used (to avoid compiler warnings).

* Added missing typecasts for some ``(void*)`` pointers (again, to avoid
  compiler warnings).

* Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of
  ``sunindextype`` in one location.

* Added missing ``#include <stdio.h>`` in ``N_Vector`` and ``SUNMatrix`` header
  files.

* Added missing prototype for :c:func:`IDASpilsGetNumJTSetupEvals`.

* Fixed an indexing bug in the CUDA ``N_Vector`` implementation of
  :c:func:`N_VWrmsNormMask` and revised the RAJA ``N_Vector`` implementation of
  :c:func:`N_VWrmsNormMask` to work with mask arrays using values other than
  zero or one. Replaced ``double`` with ``realtype`` in the RAJA vector test
  functions.

* Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not
  require a ``SUNMatrix`` module (e.g., iterative linear solvers).

In addition to the changes above, minor corrections were also made to the
example programs, build system, and user documentation.

Changes in v3.1.0
-----------------

Added ``N_Vector`` print functions that write vector data to a specified file
(e.g., :c:func:`N_VPrintFile_Serial`).

Added ``make test`` and ``make test_install`` options to the build system for
testing SUNDIALS after building with ``make`` and installing with ``make
install`` respectively.

Changes in v3.0.0
-----------------

All interfaces to matrix structures and linear solvers have been reworked, and
all example programs have been updated.  The goal of the redesign of these
interfaces was to provide more encapsulation and to ease interfacing of custom
linear solvers and interoperability with linear solver libraries.  Specific
changes include:

* Added generic ``SUNMatrix`` module with three provided implementations: dense,
  banded, and sparse. These replicate previous SUNDIALS Dls and Sls matrix
  structures in a single object-oriented API.

* Added example problems demonstrating use of generic ``SUNMatrix`` modules.

* Added generic ``SUNLinearSolver`` module with eleven provided implementations:
  SUNDIALS native dense, SUNDIALS native banded, LAPACK dense, LAPACK band, KLU,
  SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, and PCG. These replicate previous
  SUNDIALS generic linear solvers in a single object-oriented API.

* Added example problems demonstrating use of generic ``SUNLinearSolver``
  modules.

* Expanded package-provided direct linear solver (Dls) interfaces and scaled,
  preconditioned, iterative linear solver (Spils) interfaces to utilize generic
  ``SUNMatrix`` and ``SUNLinearSolver`` objects.

* Removed package-specific, linear solver-specific, solver modules
  (e.g. ``CVDENSE``, ``KINBAND``, ``IDAKLU``, ``ARKSPGMR``) since their
  functionality is entirely replicated by the generic Dls/Spils interfaces and
  ``SUNLinearSolver`` and ``SUNMatrix`` modules. The exception is ``CVDIAG``, a
  diagonal approximate Jacobian solver available to CVODE and CVODES.

* Converted all SUNDIALS example problems and files to utilize the new generic
  ``SUNMatrix`` and ``SUNLinearSolver`` objects, along with updated Dls and
  Spils linear solver interfaces.

* Added Spils interface routines to ARKODE, CVODE, CVODES, IDA, and IDAS to
  allow specification of a user-provided "JTSetup" routine.  This change
  supports users who wish to set up data structures for the user-provided
  Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup
  setup per Newton iteration can be amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added – one for CUDA and one
for RAJA vectors.  These vectors are supplied to provide very basic support for
running on GPU architectures. Users are advised that these vectors both move all
data to the GPU device upon construction, and speedup will only be realized if
the user also conducts the right-hand-side or residual function evaluation on
the device. In addition, these vectors assume the problem fits on one GPU.
For further information about RAJA, users are referred to the web site,
https://software.llnl.gov/RAJA/.  These additions are accompanied by updates
to various interface functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can
be configured to be a 32- or 64-bit integer data index type.  ``sunindextype``
is defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``.  The Fortran interfaces
continue to use ``long int`` for indices, except for their sparse matrix
interface that now uses the new ``sunindextype``.  This new flexible capability
for index types includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with
either 32-bit or 64-bit capabilities depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype``
values ``TRUE`` and ``FALSE`` have been changed to ``SUNTRUE`` and ``SUNFALSE``
respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS
type information for use in Fortran programs.

The build system was expanded to support many of the xSDK-compliant keys.  The
xSDK is a movement in scientific software to provide a foundation for the rapid
and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

Added functions :c:func:`SUNDIALSGetVersion` and
:c:func:`SUNDIALSGetVersionNumber` to get SUNDIALS release version information
at runtime.

In addition, numerous changes were made to the build system.  These include the
addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES`` CMake variables,
additional error checking during CMake configuration, minor bug fixes, and
renaming CMake options to enable/disable examples for greater clarity and an
added option to enable/disable Fortran 77 examples.  These changes included
changing ``EXAMPLES_ENABLE`` to :cmakeop:`EXAMPLES_ENABLE_C`, changing
``CXX_ENABLE`` to :cmakeop:`EXAMPLES_ENABLE_CXX`, changing ``F90_ENABLE`` to
:cmakeop:`EXAMPLES_ENABLE_F90`, and adding an :cmakeop:`EXAMPLES_ENABLE_F77`
option.

A bug fix was done to add a missing prototype for :c:func:`IDASetMaxBacksIC` in
``ida.h``.

Corrections and additions were made to the examples, to installation-related
files, and to the user documentation.

Changes in v2.9.0
-----------------

Two additional ``N_Vector`` implementations were added – one for Hypre
(parallel) ParVector vectors, and one for PETSc vectors. These additions are
accompanied by additions to various interface functions and to user
documentation.

Each ``N_Vector`` module now includes a function, :c:func:`N_VGetVectorID`, that
returns the ``N_Vector`` module name.

An optional input function was added to set a maximum number of linesearch
backtracks in the initial condition calculation.  Also, corrections were made to
three Fortran interface functions.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in solver
``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem solution.

A memory leak was fixed in the banded preconditioner interface.  In addition,
updates were done to return integers from linear solver and preconditioner
"free" functions.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant dot
product. Various additions and corrections were made to the interfaces to the
sparse solvers KLU and SuperLU_MT, including support for CSR format when using
KLU.

New examples were added for use of the OpenMP vector.

Minor corrections and additions were made to the IDA solver, to the Fortran
interfaces, to the examples, to installation-related files, and to the user
documentation.

Changes in v2.8.0
-----------------

Two major additions were made to the linear system solvers that are available
for use with the IDA solver. First, in the serial case, an interface to the
sparse direct solver KLU was added.  Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct
solver option, to be used with the serial version of the ``N_Vector`` module.
As part of these additions, a sparse matrix (CSC format) structure was added to
IDA.

Otherwise, only relatively minor modifications were made to IDA:

In :c:func:`IDARootfind`, a minor bug was corrected, where the input array
``rootdir`` was ignored, and a line was added to break out of root-search loop
if the initial interval size is below the tolerance ``ttol``.

In ``IDALapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to ``smu =
mu + ml`` to correct an illegal input error for ``DGBTRF/DGBTRS``.

A minor bug was fixed regarding the testing of the input ``tstop`` on the first
call to :c:func:`IDASolve`.

In order to avoid possible name conflicts, the mathematical macro and function
names ``MIN``, ``MAX``, ``SQR``, ``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and
``RPowerR`` were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``, ``SUNRabs``,
``SUNRsqrt``, ``SUNRexp``, ``SRpowerI``, and ``SUNRpowerR``, respectively.
These names occur in both the solver and in various example programs.

In the FIDA optional input routines ``FIDASETIIN``, ``FIDASETRIN``, and
``FIDASETVIN``, the optional fourth argument ``key_length`` was removed, with
hardcoded key string lengths passed to all ``strncmp`` tests.

In all FIDA examples, integer declarations were revised so that those which must
match a C type ``long int`` are declared ``INTEGER*8``, and a comment was added
about the type match. All other integer declarations are just
``INTEGER``. Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing
environments — one for OpenMP, denoted :ref:`NVECTOR_OPENMP <NVectors.OpenMP>`,
and one for Pthreads, denoted :ref:`NVECTOR_PTHREADS <NVectors.Pthreads>`.

With this version of SUNDIALS, support and documentation of the Autotools mode
of installation is being dropped, in favor of the CMake mode, which is
considered more widely portable.

Changes in v2.7.0
-----------------

One significant design change was made with this release: The problem size and
its relatives, bandwidth parameters, related internal indices, pivot arrays, and
the optional output ``lsflag`` have all been changed from type ``int`` to type
``long int``, except for the problem size and bandwidths in user calls to
routines specifying BLAS/LAPACK routines for the dense/band linear solvers. The
function ``NewIntArray`` is replaced by a pair ``NewIntArray`` and
``NewLintArray``, for ``int`` and ``long int`` arrays, respectively.

A large number of minor errors have been fixed. Among these are the following:
After the solver memory is created, it is set to zero before being filled.  To
be consistent with IDAS, IDA uses the function ``IDAGetDky`` for optional output
retrieval.  In each linear solver interface function, the linear solver memory
is freed on an error return, and the ``**Free`` function now includes a line
setting to NULL the main memory pointer to the linear solver memory.  A memory
leak was fixed in two of the ``IDASp***Free`` functions.  In the rootfinding
functions ``IDARcheck1`` and ``IDARcheck2``, when an exact zero is found, the
array ``glo`` of :math:`g` values at the left endpoint is adjusted, instead of
shifting the :math:`t` location ``tlo`` slightly.  In the installation files, we
modified the treatment of the macro SUNDIALS_USE_GENERIC_MATH, so that the
parameter GENERIC_MATH_LIB is either defined (with no value) or not defined.

Changes in v2.6.0
-----------------

Two new features were added in this release: (a) a new linear solver module,
based on BLAS and LAPACK for both dense and banded matrices, and (b) option to
specify which direction of zero-crossing is to be monitored while performing
rootfinding.

The user interface has been further refined. Some of the API changes involve:
(a) a reorganization of all linear solver modules into two families (besides the
already present family of scaled preconditioned iterative linear solvers, the
direct solvers, including the new LAPACK-based ones, were also organized into a
*direct* family); (b) maintaining a single pointer to user data, optionally
specified through a ``Set``-type function; (c) a general streamlining of the
band-block-diagonal preconditioner module distributed with the solver.

Changes in v2.5.0
-----------------

The main changes in this release involve a rearrangement of the entire SUNDIALS
source tree (see :numref:`IDA.Organization`). At the user interface level, the main
impact is in the mechanism of including SUNDIALS header files which must now
include the relative path (e.g. ``#include <cvode/cvode.h>``).  Additional
changes were made to the build system: all exported header files are now
installed in separate subdirectories of the installation *include* directory.

A bug was fixed in the internal difference-quotient dense and banded Jacobian
approximations, related to the estimation of the perturbation (which could have
led to a failure of the linear solver when zero components with sufficiently
small absolute tolerances were present).

The user interface to the consistent initial conditions calculations was
modified.  The :c:func:`IDACalcIC` arguments ``t0``, ``yy0``, and ``yp0`` were
removed and a new function, :c:func:`IDAGetConsistentIC` is provided.

The functions in the generic dense linear solver (``sundials_dense`` and
``sundials_smalldense``) were modified to work for rectangular :math:`m \times
n` matrices (:math:`m \le n`), while the factorization and solution functions
were renamed to ``DenseGETRF / denGETRF`` and ``DenseGETRS / denGETRS``,
respectively.  The factorization and solution functions in the generic band
linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.

Changes in v2.4.0
-----------------

FIDA, a Fortran-C interface module, was added.

IDASPBCG and IDASPTFQMR modules have been added to interface with the Scaled
Preconditioned Bi-CGstab (SPBCG) and Scaled Preconditioned Transpose-Free
Quasi-Minimal Residual (SPTFQMR) linear solver modules, respectively (for
details see :numref:IDA.Usage.CC).  At the same time, function type names for Scaled
Preconditioned Iterative Linear Solvers were added for the user-supplied
Jacobian-times-vector and preconditioner setup and solve functions.

The rootfinding feature was added, whereby the roots of a set of given functions
may be computed during the integration of the DAE system.

A user-callable routine was added to access the estimated local error vector.

The deallocation functions now take as arguments the address of the respective
memory block pointer.

To reduce the possibility of conflicts, the names of all header files have been
changed by adding unique prefixes (``ida_`` and ``sundials_``).  When using the
default installation procedure, the header files are exported under various
subdirectories of the target ``include`` directory. For more details see
Appendix :numref:`Installation`.

Changes in v2.3.0
-----------------

The user interface has been further refined. Several functions used for setting
optional inputs were combined into a single one.  An optional user-supplied
routine for setting the error weight vector was added.  Additionally, to resolve
potential variable scope issues, all SUNDIALS solvers release user data right
after its use. The build systems has been further improved to make it more
robust.

Changes in v2.2.2
-----------------

Minor corrections and improvements were made to the build system.  A new chapter
in the User Guide was added — with constants that appear in the user interface.

Changes in v2.2.1
-----------------

The changes in this minor SUNDIALS release affect only the build system.

Changes in v2.2.0
-----------------

The major changes from the previous version involve a redesign of the user
interface across the entire SUNDIALS suite. We have eliminated the mechanism of
providing optional inputs and extracting optional statistics from the solver
through the ``iopt`` and ``ropt`` arrays. Instead, IDA now provides a set of
routines (with prefix ``IDASet``) to change the default values for various
quantities controlling the solver and a set of extraction routines (with prefix
``IDAGet``) to extract statistics after return from the main solver routine.
Similarly, each linear solver module provides its own set of ``Set``- and
``Get``-type routines. For more details see :numref:`IDA.Usage.CC.optional_output`.

Additionally, the interfaces to several user-supplied routines (such as those
providing Jacobians and preconditioner information) were simplified by reducing
the number of arguments. The same information that was previously accessible
through such arguments can now be obtained through ``Get``-type functions.

Installation of IDA (and all of SUNDIALS) has been completely redesigned and is
now based on configure scripts.

.. _IDA.Introduction.Reading:

Reading this User Guide
=======================

The structure of this document is as follows:

* In Chapter :numref:`IDA.Mathematics`, we give short descriptions of the numerical
  methods implemented by IDA for the solution of initial value problems for
  systems of DAEs, along with short descriptions of preconditioning
  (:numref:`IDA.Mathematics.Preconditioning`) and rootfinding
  (:numref:`IDA.Mathematics.rootfinding`).

* The following chapter describes the structure of the SUNDIALS suite of solvers
  (:numref:`IDA.Organization`) and the software organization of the IDA solver
  (:numref:`IDA.Organization.IDA`).

* Chapter :numref:`IDA.Usage.CC` is the main usage document for IDA for C and C++
  applications. It includes a complete description of the user interface for the
  integration of DAE initial value problems. This is followed by documentation
  for using IDA with Fortran applications and on GPU accelerated systems.

* Chapter :numref:`NVectors` gives a brief overview of the generic ``N_Vector``
  module shared among the various components of SUNDIALS, as well as details on
  the ``N_Vector`` implementations provided with SUNDIALS.

* Chapter :numref:`SUNMatrix` gives a brief overview of the generic
  ``SUNMatrix`` module shared among the various components of SUNDIALS, and
  details on the ``SUNMatrix`` implementations provided with SUNDIALS.

* Chapter :numref:`SUNLinSol` gives a brief overview of the generic
  ``SUNLinearSolver`` module shared among the various components of
  SUNDIALS. This chapter contains details on the ``SUNLinearSolver``
  implementations provided with SUNDIALS.  The chapter also contains details on
  the ``SUNLinearSolver`` implementations provided with SUNDIALS that interface
  with external linear solver libraries.

* Chapter :numref:`SUNNonlinSol` describes the ``SUNNonlinearSolver`` API and
  nonlinear solver implementations shared among the various components of
  SUNDIALS.

* Finally, in the appendices, we provide detailed instructions for the
  installation of IDA, within the structure of SUNDIALS (Appendix
  :numref:`Installation`), as well as a list of all the constants used for input
  to and output from IDA functions (Appendix :numref:`IDA.Constants`).

..
   Finally, the reader should be aware of the following notational conventions in
   this user guide: program listings and identifiers (such as :c:func:`IDAInit`)
   within textual explanations appear in typewriter type style; fields in C
   structures (such as *content*) appear in italics; and packages or modules, such
   as IDADLS, are written in all capitals.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.
