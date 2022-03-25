.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODE.Introduction:

************
Introduction
************

CVODE is part of a software family called SUNDIALS: SUite of
Nonlinear and DIfferential/ALgebraic equation
Solvers :cite:p:`HBGLSSW:05`. This suite consists of CVODE,
ARKODE, KINSOL, and IDA, and variants of these with sensitivity
analysis capabilities.

.. _CVODE.Introduction.history:

Historical Background
=====================

Fortran solvers for ODE initial value problems are widespread and heavily
used. Two solvers that have been written at LLNL in the past are
VODE :cite:p:`BBH:89` and
VODPK :cite:p:`Byr:92`. VODE is a general purpose
solver that includes methods for both stiff and nonstiff systems, and in
the stiff case uses direct methods (full or banded) for the solution of
the linear systems that arise at each implicit step. Externally,
VODE is very similar to the well known solver
LSODE :cite:p:`RaHi:94`. VODPK is a variant of VODE
that uses a preconditioned Krylov (iterative) method, namely GMRES, for
the solution of the linear systems. VODPK is a powerful tool for
large stiff systems because it combines established methods for stiff
integration, nonlinear iteration, and Krylov (linear) iteration with a
problem-specific treatment of the dominant source of stiffness, in the
form of the user-supplied preconditioner
matrix :cite:p:`BrHi:89`. The capabilities of both VODE
and VODPK have been combined in the C-language package
CVODE :cite:p:`CoHi:96`.

At present, CVODE may utilize a variety of Krylov methods provided
in SUNDIALS that can be used in conjuction with Newton iteration:
these include the GMRES (Generalized Minimal
RESidual) :cite:p:`SaSc:86`, FGMRES (Flexible Generalized
Minimum RESidual) :cite:p:`Saa:93`, Bi-CGStab (Bi-Conjugate
Gradient Stabilized) :cite:p:`Van:92`, TFQMR (Transpose-Free
Quasi-Minimal Residual) :cite:p:`Fre:93`, and PCG
(Preconditioned Conjugate Gradient) :cite:p:`HeSt:52` linear
iterative methods. As Krylov methods, these require almost no matrix
storage for solving the Newton equations as compared to direct methods.
However, the algorithms allow for a user-supplied preconditioner matrix,
and for most problems preconditioning is essential for an efficient
solution. For very large stiff ODE systems, the Krylov methods are
preferable over direct linear solver methods, and are often the only
feasible choice. Among the Krylov methods in SUNDIALS, we recommend
GMRES as the best overall choice. However, users are encouraged to
compare all options, especially if encountering convergence failures
with GMRES. Bi-CGStab and TFQMR have an advantage in storage
requirements, in that the number of workspace vectors they require is
fixed, while that number for GMRES depends on the desired Krylov
subspace size. FGMRES has an advantage in that it is designed to support
preconditioners that vary between iterations (e.g. iterative methods).
PCG exhibits rapid convergence and minimal workspace vectors, but only
works for symmetric linear systems.

In the process of translating the VODE and VODPK algorithms into
C, the overall CVODE organization has been changed considerably.
One key feature of the CVODE organization is that the linear system
solvers comprise a layer of code modules that is separated from the
integration algorithm, allowing for easy modification and expansion of
the linear solver array. A second key feature is a separate module
devoted to vector operations; this facilitated the extension to
multiprosessor environments with minimal impacts on the rest of the
solver, resulting in PVODE :cite:p:`ByHi:99`, the parallel
variant of CVODE.

Around 2002, the functionality of CVODE and PVODE were combined
into one single code, simply called CVODE. Development of this
version of CVODE was concurrent with a redesign of the vector
operations module across the SUNDIALS suite. The key feature of the
``N_Vector`` module is that it is written in terms of abstract vector
operations with the actual vector kernels attached by a particular
implementation (such as serial or parallel) of ``N_Vector``. This allows
writing the SUNDIALS solvers in a manner independent of the actual
``N_Vector`` implementation (which can be user-supplied), as well as
allowing more than one ``N_Vector`` module linked into an executable
file. SUNDIALS (and thus CVODE) is supplied with a wide range of different
``N_Vector`` implementations, including: serial, MPI-parallel, both OpenMP and
Pthreads thread-parallel ``N_Vector`` implementations, a Hypre parallel
implementation, a PETSc implementation, and various GPU-enabled
implementations.

.. There are several motivations for choosing the C language for
.. CVODE. First, a general movement away from Fortran and toward C in
.. scientific computing was apparent. Second, the pointer, structure, and
.. dynamic memory allocation features in C are extremely useful in software
.. of this complexity, with the great variety of method options offered.
.. Finally, we prefer C over |CPP| for CVODE because of the
.. wider availability of C compilers, the potentially greater
.. efficiency of C, and the greater ease of interfacing the solver to
.. applications written in extended Fortran.

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

**CVODE**

The previously deprecated function ``CVodeSetMaxStepsBetweenJac`` has been
removed and replaced with :c:func:`CVodeSetJacEvalFrequency`.

The CVODE Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
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
| ``CVSpilsSetLinearSolver``      | ``CVodeSetLinearSolver``       |
+---------------------------------+--------------------------------+
| ``CVSpilsSetEpsLin``            | ``CVodeSetEpsLin``             |
+---------------------------------+--------------------------------+
| ``CVSpilsSetPreconditioner``    | ``CVodeSetPreconditioner``     |
+---------------------------------+--------------------------------+
| ``CVSpilsSetJacTimes``          | ``CVodeSetJacTimes``           |
+---------------------------------+--------------------------------+
| ``CVSpilsGetWorkSpace``         | ``CVodeGetLinWorkSpace``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumPrecEvals``      | ``CVodeGetNumPrecEvals``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumPrecSolves``     | ``CVodeGetNumPrecSolves``      |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumLinIters``       | ``CVodeGetNumLinIters``        |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumConvFails``      | ``CVodeGetNumConvFails``       |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumJTSetupEvals``   | ``CVodeGetNumJTSetupEvals``    |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumJtimesEvals``    | ``CVodeGetNumJtimesEvals``     |
+---------------------------------+--------------------------------+
| ``CVSpilsGetNumRhsEvals``       | ``CVodeGetNumLinRhsEvals``     |
+---------------------------------+--------------------------------+
| ``CVSpilsGetLastFlag``          | ``CVodeGetLastLinFlag``        |
+---------------------------------+--------------------------------+
| ``CVSpilsGetReturnFlagName``    | ``CVodeGetLinReturnFlagName``  |
+---------------------------------+--------------------------------+
| ``CVDlsSetLinearSolver``        | ``CVodeSetLinearSolver``       |
+---------------------------------+--------------------------------+
| ``CVDlsSetJacFn``               | ``CVodeSetJacFn``              |
+---------------------------------+--------------------------------+
| ``CVDlsGetWorkSpace``           | ``CVodeGetLinWorkSpace``       |
+---------------------------------+--------------------------------+
| ``CVDlsGetNumJacEvals``         | ``CVodeGetNumJacEvals``        |
+---------------------------------+--------------------------------+
| ``CVDlsGetNumRhsEvals``         | ``CVodeGetNumLinRhsEvals``     |
+---------------------------------+--------------------------------+
| ``CVDlsGetLastFlag``            | ``CVodeGetLastLinFlag``        |
+---------------------------------+--------------------------------+
| ``CVDlsGetReturnFlagName``      | ``CVodeGetLinReturnFlagName``  |
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
support the SYCL backend in addition to the CUDA and HIP backend. Users can
choose the backend when configuring SUNDIALS by using the
``SUNDIALS_RAJA_BACKENDS`` CMake variable.  This module remains experimental
and is subject to change from version to version.

New :c:type:`SUNMatrix` and :c:type:`SUNLinearSolver` implementations were added to
interface with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix
and the linear solver support general dense linear systems as well as block
diagonal linear systems. See :numref:`SUNLinSol.OneMklDense` for more details.
This module is experimental and is subject to change from version to version.

Added a new *optional* function to the SUNLinearSolver API,
:c:func:`SUNLinSolSetZeroGuess()`, to indicate that the next call to
:c:func:`SUNlinSolSolve()` will be made with a zero initial guess. SUNLinearSolver
implementations that do not use the :c:func:`SUNLinSolNewEmpty` constructor will,
at a minimum, need set the ``setzeroguess`` function pointer in the linear solver
``ops`` structure to ``NULL``. The SUNDIALS iterative linear solver
implementations have been updated to leverage this new set function to remove
one dot product per solve.

CVODE now supports a new "matrix-embedded" :c:type:`SUNLinearSolver` type.  This
type supports user-supplied :c:type:`SUNLinearSolver` implementations that set up
and solve the specified linear system at each linear solve call.  Any
matrix-related data structures are held internally to the linear solver itself,
and are not provided by the SUNDIALS package.

Added specialized fused HIP kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the
:ref:`N_Vector HIP <NVectors.HIP>` module. See the optional input function
:c:func:`CVodeSetUseIntegratorFusedKernels()` for more information. As with
other SUNDIALS HIP features, this capability is considered experimental and may
change from version to version.

Added the function :c:func:`CVodeSetNlsRhsFn()` to supply an alternative right-hand
side function for use within nonlinear system function evaluations.

The installed ``SUNDIALSConfig.cmake`` file now supports the ``COMPONENTS`` option
to ``find_package``. The exported targets no longer have ``IMPORTED_GLOBAL``
set.

A bug was fixed in :c:func:`SUNMatCopyOps()` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the :c:ref:`SPBCGS <SUNLinSol.SPBCGS>` and
:c:ref:`SPTFQMR <SUNLinSol.SPTFQMR>` solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using :c:ref:`SPBCGS <SUNLinSol.SPBCGS>` or :c:ref:`SPTFQMR <SUNLinSol.SPTFQMR>`
as standalone solvers as all SUNDIALS packages utilize a zero initial guess.

Changes in v5.7.0
-----------------

A new :c:type:`N_Vector` implementation based on the SYCL abstraction layer
has been added targeting Intel GPUs. At present the only SYCL
compiler supported is the DPC++ (Intel oneAPI) compiler. See
:numref:`NVectors.sycl` for more details. This module is
considered experimental and is subject to major changes even in minor
releases.

New ``SUNMatrix`` and ``SUNLinearSolver`` implementations were added to
interface with the MAGMA linear algebra library. Both the matrix and the
linear solver support general dense linear systems as well as block
diagonal linear systems, and both are targeted at GPUs (AMD or NVIDIA).
See :numref:`SUNLinSol.MagmaDense` for more details.

Changes in v5.6.1
-----------------

Fixed a bug in the SUNDIALS CMake which caused an error if the
``CMAKE_CXX_STANDARD`` and ``SUNDIALS_RAJA_BACKENDS`` options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

Changes in v5.6.0
-----------------

A new :c:type:`N_Vector` implementation based on the AMD ROCm HIP platform has
been added. This vector can target NVIDIA or AMD GPUs. See
:numref:`NVectors.hip` for more details. This module is
considered experimental and is subject to change from version to
version.

The :ref:`RAJA N_Vector <NVectors.RAJA>` implementation has been updated to support the HIP
backend in addition to the CUDA backend. Users can choose the backend
when configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the N_Vector API. This
operation is useful for N_Vectors that utilize dual memory spaces, e.g.
the native SUNDIALS CUDA N_Vector.

The :ref:`SUNMATRIX_CUSPARSE <SUNMatrix.cuSparse>` and
:ref:`SUNLINEARSOLVER_CUSOLVERSP_BATCHQR <SUNLinSol.cuSolverSp>`
implementations no longer require the SUNDIALS :ref:`CUDA N_Vector <NVectors.CUDA>`. Instead,
they require that the vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and that
the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA device pointer.

Changes in v5.5.0
-----------------

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now
required. Users will likely see deprecation warnings, but otherwise the
changes should be fully backwards compatible for almost all users.
SUNDIALS now exports CMake targets and installs a
SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.

Changes in v5.4.0
-----------------

Added new functions :c:func:`CVodeComputeState`, and
:c:func:`CVodeGetNonlinearSystemData` which advanced users might find useful if
providing a custom :c:type:`SUNNonlinSolSysFn`.

Added the function :c:func:`CVodeSetLSNormFactor` to specify the factor for
converting between integrator tolerances (WRMS norm) and linear solver
tolerances (L2 norm) i.e., ``tol_L2 = nrmfac * tol_WRMS``.

The expected behavior of :c:func:`SUNNonlinSolGetNumIters` and
:c:func:`SUNNonlinSolGetNumConvFails` in the :c:type:`SUNNonlinearSolver` API have
been updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided :c:type:`SUNNonlinearSolver`
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retreived by calling
:c:func:`CVodeGetNumNonlinSolvIters`, the cumulative number of failures with
:c:func:`CVodeGetNumNonlinSolvConvFails`, or both with
:c:func:`CVodeGetNonlinSolvStats`.

A minor inconsistency in checking the Jacobian evaluation frequency has
been fixed. As a result codes using using a non-default Jacobian update
frequency through a call to :c:func:`CVodeSetMaxStepsBetweenJac` will need to increase the provided value by
1 to achieve the same behavior as before. For greater clarity the
function has been deprecated and replaced with :c:func:`CVodeSetJacEvalFrequency`. Additionally, the
function :c:func:`CVodeSetLSetupFrequency` has been added to set the frequency of calls to the linear
solver setup function.

A new class, :ref:`SUNMemoryHelper <SUNMemory>`, was added to support **GPU
users** who have complex memory management needs such as using memory pools.
This is paired with new constructors for the ``NVECTOR_CUDA`` and
``NVECTOR_RAJA`` modules that accept a ``SUNMemoryHelper`` object. Refer to
:numref:`SUNDIALS.GPU`, :numref:`SUNMemory`, :numref:`NVectors.cuda` and
:numref:`NVectors.raja` for more information.

The ``NVECTOR_RAJA`` vector implementation has been updated to mirror the
``NVECTOR_CUDA`` implementation. Notably, the update adds managed memory
support. Users of the vector will need to update any calls to the function
because that signature was changed. This vector remains experimental and is
subject to change from version to version.

The ``NVECTOR_TRILINOS`` vector implementation has been updated to work with
Trilinos 12.18+. This update changes the local ordinal type to always be an
``int``.

Changes in v5.3.0
-----------------

Fixed a bug in the iterative linear solver modules where an error is not
returned if the Atimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added specialized fused CUDA kernels to CVODE which may offer
better performance on smaller problems when using CVODE with the
``NVECTOR_CUDA`` module. See the optional input function for more
information. As with other SUNDIALS CUDA features, this
capability is considered experimental and may change from version to
version.

Added the ability to control the CUDA kernel launch parameters for the
``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE`` modules. These modules remain
experimental and are subject to change from version to version. In addition, the
kernels were rewritten to be more flexible. Most users should see equivalent
performance or some improvement, but a select few may observe minor performance
degradation with the default settings. Users are encouraged to contact the
SUNDIALS team about any perfomance changes that they notice.

Added new capabilities for monitoring the solve phase in the
``SUNNONLINSOL_NEWTON`` and ``SUNNONLINSOL_FIXEDPOINT`` modules, and the
SUNDIALS iterative linear solver modules. SUNDIALS must be built
with the ``SUNDIALS_BUILD_WITH_MONITORING`` CMake option set to ``TRUE`` to use these capabilties.

Added a new function, :c:func:`CVodeSetMonitorFn`, that takes a user-function to be called by
CVODE after every :math:`nst` succesfully completed time-steps. This
is intended to provide a way of monitoring the CVODE statistics
throughout the simulation.

Added a new function :c:func:`CVodeGetLinSolveStats` to get the CVODE linear solver statistics as a
group.

Added the optional function :c:func:`CVodeSetJacTimsRhsFn` to specify an alternative right-hand side
function for computing Jacobian-vector products with the internal
difference quotient approximation.

Added support for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user
defined projection function. This implementation is accompanied by
additions to user documentation and CVODE examples. See
:c:func:`CVodeSetConstraints` for more information.

Added support for CUDA v11.

Changes in v5.2.0
-----------------

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``, ``xlf2003``,
or ``xlf2003_r``.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attributes missing on some SUNDIALS API functions.

Added a new ``SUNMatrix`` implementation, ``SUNMATRIX_CUSPARSE``, that
interfaces to the sparse matrix implementation from the NVIDIA cuSPARSE library.
In addition, the linear solver has been updated to use this matrix, therefore,
users of this module will need to update their code. These modules are still
considered to be experimental, thus they are subject to breaking changes even in
minor releases.

The function :c:func:`CVodeSetLinearSolutionScaling` was added to enable or
disable the scaling applied to linear system solutions with matrix-based linear
solvers to account for a lagged value of :math:`\gamma` in the linear system
matrix :math:`I - \gamma J`. Scaling is enabled by default when using a
matrix-based linear solver with BDF methods.

Changes in v5.1.0
-----------------

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables ``PETSC_INCLUDES`` and ``PETSC_LIBRARIES`` instead of ``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify the CUDA
architecture to compile for.

Added two utility functions, :c:func:`SUNDIALSFileOpen` and
:c:func:`SUNDIALSFileClose` for creating/destroying file pointers that are
useful when using the Fortran 2003 interfaces.

Added support for constant damping to the :ref:`SUNNonlinearSolver_FixedPoint
<SUNNonlinSol.FixedPoint>` module when using Anderson acceleration.

Changes in v5.0.0
-----------------

**Build system changes**

-  Increased the minimum required CMake version to 3.5 for most
   SUNDIALS configurations, and 3.10 when CUDA or OpenMP with device
   offloading are enabled.

-  The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been removed to simplify
   builds as SUNDIALS packages do not use BLAS directly. For third
   party libraries that require linking to BLAS, the path to the BLAS
   library should be included in the variable for the third party
   library *e.g.*, ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

-  Fixed a bug in the build system that prevented the ``NVECTOR_PTHREADS``
   module from being built.

**NVECTOR module changes**

-  Two new functions were added to aid in creating custom ``N_Vector``
   objects. The constructor :c:func:`N_VNewEmpty` allocates an “empty” generic ``N_Vector``
   with the object’s content pointer and the function pointers in the
   operations structure initialized to  ``NULL``. When used in the constructor
   for custom objects this function will ease the introduction of any
   new optional operations to the ``N_Vector`` API by ensuring only
   required operations need to be set. Additionally, the function :c:func:`N_VCopyOps` has
   been added to copy the operation function pointers between vector
   objects. When used in clone routines for custom vector objects these
   functions also will ease the introduction of any new optional
   operations to the ``N_Vector`` API by ensuring all operations are
   copied when cloning objects. See :numref:`NVectors.Description.custom_implementation` for more details.

-  Two new ``N_Vector`` implementations, ``NVECTOR_MANYVECTOR`` and
   ``NVECTOR_MPIMANYVECTOR``, have been created to support flexible
   partitioning of solution data among different processing elements
   (e.g., CPU + GPU) or for multi-physics problems that couple distinct
   MPI-based simulations together. This implementation is accompanied by
   additions to user documentation and SUNDIALS examples. See
   :numref:`NVectors.manyvector` and :numref:`NVectors.mpimanyvector` for more
   details.

-  One new required vector operation and ten new optional vector
   operations have been added to the ``N_Vector`` API. The new required
   operation, , returns the global length of an . The optional operations have
   been added to support the new ``NVECTOR_MPIMANYVECTOR`` implementation. The
   operation must be implemented by subvectors that are combined to create an
   ``NVECTOR_MPIMANYVECTOR``, but is not used outside of this context. The
   remaining nine operations are optional local reduction operations intended to
   eliminate unnecessary latency when performing vector reduction operations
   (norms, etc.) on distributed memory systems. The optional local reduction
   vector operations are :c:func:`N_VDotProdLocal`, :c:func:`N_VMaxNormLocal`,
   :c:func:`N_VL1NormLocal`, :c:func:`N_VWSqrSumLocal`,
   :c:func:`N_VWSqrSumMaskLocal`, :c:func:`N_VInvTestLocal`,
   :c:func:`N_VConstrMaskLocal`, :c:func:`N_VMinLocal`, and
   :c:func:`N_VMinQuotientLocal`. If an ``N_Vector`` implementation defines any
   of the local operations as , then the ``NVECTOR_MPIMANYVECTOR`` will call
   standard ``N_Vector`` operations to complete the computation.

-  An additional ``N_Vector`` implementation, ``NVECTOR_MPIPLUSX``, has been
   created to support the MPI+X paradigm where X is a type of on-node
   parallelism (*e.g.*, OpenMP, CUDA). The implementation is accompanied
   by additions to user documentation and SUNDIALS examples. See
   :numref:`NVectors.mpiplusx` for more details.

-  The and functions have been removed from the ``NVECTOR_CUDA`` and
   ``NVECTOR_RAJA`` implementations respectively. Accordingly, the
   ``nvector_mpicuda.h``, ``libsundials_nvecmpicuda.lib``,
   ``libsundials_nvecmpicudaraja.lib``, and files have been removed. Users
   should use the ``NVECTOR_MPIPLUSX`` module coupled in conjunction with the
   ``NVECTOR_CUDA`` or ``NVECTOR_RAJA`` modules to replace the functionality.
   The necessary changes are minimal and should require few code modifications.
   See the programs in and for examples of how to use the ``NVECTOR_MPIPLUSX``
   module with the ``NVECTOR_CUDA`` and ``NVECTOR_RAJA`` modules respectively.

-  Fixed a memory leak in the ``NVECTOR_PETSC`` module clone function.

-  Made performance improvements to the ``NVECTOR_CUDA`` module. Users who
   utilize a non-default stream should no longer see default stream
   synchronizations after memory transfers.

-  Added a new constructor to the ``NVECTOR_CUDA`` module that allows a user
   to provide custom allocate and free functions for the vector data
   array and internal reduction buffer. See :numref:`NVectors.Cuda` for more details.

-  Added new Fortran 2003 interfaces for most ``N_Vector`` modules. See
   :numref:`NVectors` for more details on how to use
   the interfaces.

-  Added three new ``N_Vector`` utility functions :c:func:`N_VGetVecAtIndexVectorArray`,
   :c:func:`N_VSetVecAtIndexVectorArray`, and :c:func:`N_VNewVectorArray` for working
   with arrays when using the Fortran 2003 interfaces.

**SUNMatrix module changes**

-  Two new functions were added to aid in creating custom ``SUNMatrix``
   objects. The constructor :c:func:`SUNMatNewEmpty` allocates an “empty” generic ``SUNMatrix``
   with the object’s content pointer and the function pointers in the
   operations structure initialized to . When used in the constructor
   for custom objects this function will ease the introduction of any
   new optional operations to the ``SUNMatrix`` API by ensuring only
   required operations need to be set. Additionally, the function :c:func:`SUNMatCopyOps` has
   been added to copy the operation function pointers between matrix
   objects. When used in clone routines for custom matrix objects these
   functions also will ease the introduction of any new optional
   operations to the ``SUNMatrix`` API by ensuring all operations are
   copied when cloning objects. See :numref:`SUNMatrix` for more
   details.
-  A new operation, :c:func:`SUNMatMatvecSetup`, was added to the ``SUNMatrix`` API to perform any
   setup necessary for computing a matrix-vector product. This operation
   is useful for ``SUNMatrix`` implementations which need to prepare the
   matrix itself, or communication structures before performing the
   matrix-vector product. Users who have implemented custom
   ``SUNMatrix`` modules will need to at least update their code to set
   the corresponding structure member to ``NULL``. See :numref:`SUNMatrix.Ops`
   for more details.
-  The generic ``SUNMatrix`` API now defines error codes to be returned
   by ``SUNMatrix`` operations. Operations which return an integer flag
   indiciating success/failure may return different values than
   previously. See :numref:`SUNMatrix.Ops.errorCodes` for
   more details.
-  A new ``SUNMatrix`` (and ``SUNLinearSolver``) implementation was added to
   facilitate the use of the SuperLU_DIST library with SUNDIALS. See
   :numref:`SUNMatrix.SLUNRloc` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNMatrix`` modules. See
   :numref:`SUNMatrix` for more details on how to
   use the interfaces.

**SUNLinearSolver module changes**

-  A new function was added to aid in creating custom ``SUNLinearSolver``
   objects. The constructor allocates an “empty” generic ``SUNLinearSolver``
   with the object’s content pointer and the function pointers in the operations
   structure initialized to . When used in the constructor for custom objects
   this function will ease the introduction of any new optional operations to
   the ``SUNLinearSolver`` API by ensuring only required operations need to be
   set. See :numref:`SUNLinSol.API.Custom` for more details.
-  The return type of the ``SUNLinearSolver`` API function has changed from to
   to be consistent with the type used to store row indices in dense and banded
   linear solver modules.
-  Added a new optional operation to the ``SUNLinearSolver`` API,
   :c:func:`SUNLinSolLastFlag`, that returns a for identifying the linear solver module.
-  The ``SUNLinearSolver`` API has been updated to make the initialize and
   setup functions optional.
-  A new ``SUNLinearSolver`` (and ``SUNMatrix``) implementation was added to
   facilitate the use of the SuperLU_DIST library with SUNDIALS. See
   :numref:`SUNLinSol.SuperLUDIST` for more details.
-  Added a new ``SUNLinearSolver`` implementation, :ref:`SUNLINEARSOLVER_CUSOLVERSP <SUNLinSol.cuSolverSp>`,
   which leverages the NVIDIA cuSOLVER sparse batched QR method for efficiently solving block
   diagonal linear systems on NVIDIA GPUs.
-  Added three new accessor functions to the ``SUNLINSOL_KLU`` module, :c:func:`SUNLinSol_KLUGetSymbolic`,
   , :c:func:`SUNLinSol_KLUGetNumeric` and :c:func:`SUNLinSol_KLUGetCommon`, to
   provide user access to the underlying KLU solver structures. See
   :numref:`SUNLinSol.KLU` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNLinearSolver`` modules. See
   :numref:`SUNLinSol` for more details on how to use the interfaces.

**SUNNonlinearSolver module changes**

-  A new function was added to aid in creating custom ``SUNNonlinearSolver``
   objects. The constructor :c:func:`SUNNonlinSolSetConvTestFN` allocates an
   “empty” generic ``SUNNonlinearSolver`` with the object’s content pointer and
   the function pointers in the operations structure initialized to . When used
   in the constructor for custom objects this function will ease the
   introduction of any new optional operations to the ``SUNNonlinearSolver`` API
   by ensuring only required operations need to be set. See
   :numref:`SUNNonlinSol.API.Custom` for more details.
-  To facilitate the use of user supplied nonlinear solver convergence
   test functions the function in the ``SUNNonlinearSolver`` API has been
   updated to take a data pointer as input. The supplied data pointer will be
   passed to the nonlinear solver convergence test function on each call.
-  The inputs values passed to the first two inputs of the function
   :c:func:`SUNNonlinSolSolve` in the ``SUNNonlinearSolver`` have been changed to
   be the predicted state and the initial guess for the correction to that state. Additionally, the
   definitions of :c:func:`SUNNonlinSolLSetupFn` and
   :c:func:`SUNNonlinSolLSolveFn` in the ``SUNNonlinearSolver`` API have been
   updated to remove unused input parameters. For more information on the
   nonlinear system formulation see :numref:`SUNNonlinSol.CVODE` and for more
   details on the API functions see :numref:`SUNNonlinSol`.
-  Added a new ``SUNNonlinearSolver`` implementation, ``SUNNONLINSOL_PETSC``,
   which interfaces to the PETSc SNES nonlinear solver API. See
   :numref:`SUNNonlinSol.PetscSNES` for more details.
-  Added new Fortran 2003 interfaces for most ``SUNNonlinearSolver`` modules.
   See :numref:`SUNDIALS.Fortran` for more details on how to use the
   interfaces.

**CVODE changes**

-  Fixed a bug in the CVODE constraint handling where the step size
   could be set below the minimum step size.
-  Fixed a bug in the CVODE nonlinear solver interface where the
   norm of the accumulated correction was not updated when using a non-default
   convergence test function.
-  Fixed a memeory leak in FCVODE when not using the default
   nonlinear solver.
-  Removed extraneous calls to for simulations where the scalar valued
   absolute tolerance, or all entries of the vector-valued absolute tolerance
   array, are strictly positive. In this scenario, CVODE will remove at least
   one global reduction per time step.
-  The CVLS interface has been updated to only zero the Jacobian matrix
   before calling a user-supplied Jacobian evaluation function when the attached
   linear solver has type ``SUNLINEARSOLVER_DIRECT``.
-  A new linear solver interface function :c:func:`CVLsLinSysFn` was added as an
   alternative method for evaluating the linear system :math:`M = I - \gamma J`.
-  Added two new functions, :c:func:`CVodeGetCurrentGamma` and :c:func:`CVodeGetCurrentState`, which may be useful to users who
   choose to provide their own nonlinear solver implementations.
-  The CVODE Fortran 2003 interface was completely redone to be more
   sustainable and to allow users to write more idiomatic Fortran. See
   :numref:`SUNDIALS.Fortran` for more details.

Changes in v4.1.0
-----------------

An additional ``N_Vector`` implementation was added for the Tpetra
vector from the Trilinos library to facilitate interoperability
between SUNDIALS and Trilinos. This implementation is accompanied
by additions to user documentation and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in
some use cases.

The CMake option ``EXAMPLES_ENABLE_RAJA`` has been removed. The option enables all examples that
use CUDA including the RAJA examples with a CUDA back end (if the RAJA
``N_Vector`` is enabled).

The implementation header file is no longer installed. This means users
who are directly manipulating the structure will need to update their
code to use CVODE’s public API.

Python is no longer required to run ``make test`` and ``make test_install``.

Changes in v4.0.2
-----------------

Added information on how to contribute to SUNDIALS and a
contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to
a source file. The symbols are now included in the CVODE library, ``libsundials_cvode``.

Changes in v4.0.1
-----------------

No changes were made in this release.

Changes in v4.0.0
-----------------

CVODE’s previous direct and iterative linear solver interfaces,
CVDLS and CVSPILS, have been merged into a single unified linear
solver interface, CVLS, to support any valid ``SUNLinearSolver`` module.
This includes the “DIRECT” and “ITERATIVE” types as well as the new
“MATRIX_ITERATIVE” type. Details regarding how CVLS utilizes linear
solvers of each type as well as discussion regarding intended use cases
for user-supplied ``SUNLinearSolver`` implementations are included in
:numref:`SUNLinSol`. All CVODE example programs
and the standalone linear solver examples have been updated to use the
unified linear solver interface.

The unified interface for the new CVLS module is very similar to the
previous CVDLS and CVSPILS interfaces. To minimize challenges in
user migration to the new names, the previous C and Fortran routine
names may still be used; these will be deprecated in future releases, so
we recommend that users migrate to the new names soon. Additionally, we
note that Fortran users, however, may need to enlarge their array of
optional integer outputs, and update the indices that they query for
certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided ``SUNLinearSolver``
implementations have been updated to follow the naming convention
``SUNLinSol_*`` where is the name of the linear solver. Solver-specific “set”
routine names have been similarly standardized. To minimize challenges in user
migration to the new names, the previous routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. All CVODE example programs and the standalone linear solver
examples have been updated to use the new naming convention.

The :ref:`SUNMATRIX_BAND <SUNMatrix.Band>` constructor has been simplified to remove the storage upper
bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear
solver modules defined through the ``SUNNonlinearSolver`` API. This API will
ease the addition of new nonlinear solver options and allow for external
or user-supplied nonlinear solvers. The ``SUNNonlinearSolver`` API and
SUNDIALS provided modules are described in
:numref:`SUNNonlinSol` and follow the same
object oriented design and implementation used by the ``N_Vector``,
``SUNMatrix``, and ``SUNLinearSolver`` modules. Currently two ``SUNNonlinearSolver``
implementations are provided, ``SUNNONLINSOL_NEWTON`` and
``SUNNONLINSOL_FIXEDPOINT``. These replicate the previous integrator
specific implementations of a Newton iteration and a fixed-point
iteration (previously referred to as a functional iteration),
respectively. Note the ``SUNNONLINSOL_FIXEDPOINT`` module can optionally
utilize Anderson’s method to accelerate convergence. Example programs
using each of these nonlinear solver modules in a standalone manner have
been added and all CVODE example programs have been updated to use
generic ``SUNNonlinearSolver`` modules.

With the introduction of ``SUNNonlinearSolver`` modules, the ``iter`` input
parameter to :c:func:`CVodeCreate` has been removed along with the function
:c:func:`CVodeSetIterType` and the constants ``CV_NEWTON`` and
``CV_FUNCTIONAL``. Similarly, the parameter has been removed from the Fortran
interface function ``FCVMALLOC``. Instead of specifying the nonlinear iteration
type when creating the CVODE memory structure, CVODE uses the
``SUNNONLINSOL_NEWTON`` module implementation of a Newton iteration by default.
For details on using a non-default or user-supplied nonlinear solver see
:numref:CVODE.Usage.CC. CVODE functions for setting the nonlinear solver options
(e.g., :c:func:`CVodeSetMaxNonlinIters`) or getting nonlinear solver statistics
(e.g., :c:func:`CVodeGetNumNonlinSolvIters`) remain unchanged and internally
call generic ``SUNNonlinearSolver`` functions as needed.

Three fused vector operations and seven vector array operations have been added
to the ``N_Vector`` API. These *optional* operations are disabled by default and
may be activated by calling vector specific routines after creating an
``N_Vector`` (see :numref:`NVectors` for more details). The new operations
are intended to increase data reuse in vector operations, reduce parallel
communication on distributed memory systems, and lower the number of kernel
launches on systems with accelerators. The fused operations are
:c:func:`N_VLinearCombination`,  :c:func:`N_VScaleAddMulti`, and
:c:func:`N_VDotProdMulti` and the vector array operations are
:c:func:`N_VLinearCombinationVectorArray`, :c:func:`N_VScaleVectorArray`,
:c:func:`N_VConstVectorArray`, :c:func:`N_VWrmsNormVectorArray`,
:c:func:`N_VWrmsNormMaskVectorArray`, and :c:func:`N_VScaleAddMultiVectorArray`.
If an ``N_Vector`` implementation defines any of these operations as, then
standard ``N_Vector`` operations will automatically be called as necessary to
complete the computation.

Multiple updates to ``NVECTOR_CUDA`` were made:

* Changed to return the global vector length instead of the local
  vector length.
* Added to return the local vector length.
* Added to return the MPI communicator used.
* Removed the accessor functions in the namespace suncudavec.
* Changed the function to take a host data pointer and a device data
  pointer instead of an object.
* Added the ability to set the used for execution of the ``NVECTOR_CUDA``
  kernels. See the function :c:func:`N_VSetCudaStream_Cuda()`.
* Added :c:func:`N_VNewManaged_Cuda`, :c:func:`N_VMakeManaged_Cuda`, and
  :c:func:`N_VIsManagedMemory_Cuda()` functions to accommodate using managed
  memory with ``NVECTOR_CUDA``.

Multiple changes to ``NVECTOR_RAJA`` were made:

   - Changed to return the global vector length instead of the local vector length.
   - Added to return the local vector length.
   - Added to return the MPI communicator used.
   - Removed the accessor functions in the namespace suncudavec.
   - A new ``N_Vector`` implementation for leveraging OpenMP 4.5+ device
     offloading has been added, ``NVECTOR_OPENMPDEV``.
   - Two changes were made in the CVODE/CVODES/ARKODE initial step size algorithm:

     - Fixed an efficiency bug where an extra call to the right hand side function was made.
     - Changed the behavior of the algorithm if the max-iterations case is hit. Before the algorithm would exit with the step size calculated on the penultimate iteration. Now it will exit with the step size calculated on the final iteration.

A Fortran 2003 interface to CVODE has been added along with Fortran 2003
interfaces to the following shared SUNDIALS modules:

   -  ``SUNNONLINSOL_FIXEDPOINT`` and ``SUNNONLINSOL_NEWTON`` nonlinear solver modules
   -  ``SUNLINSOL_BAND``, ``SUNLINSOL_DENSE``, ``SUNLINSOL_KLU``, ``SUNLINSOL_PCG``, ``SUNLINSOL_SPBCGS``, ``SUNLINSOL_SPFGMR``, ``SUNLINSOL_SPGMR``, and ``SUNLINSOL_SPTFQMR`` linear solver modules
   -  ``NVECTOR_SERIAL``, ``NVECTOR_PTHREADS``, and ``NVECTOR_OPENMP`` vector modules

Changes in v3.2.1
-----------------

The changes in this minor release include the following:

-  Fixed a bug in the CUDA ``N_Vector`` where the operation could
   write beyond the allocated vector data.

-  Fixed library installation path for multiarch systems. This fix
   changes the default library installation path to
   ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`` from
   ``CMAKE_INSTALL_PREFIX/lib``. ``CMAKE_INSTALL_LIBDIR`` is automatically set,
   but is available as a CMake option that can modified.

Changes in v3.2.0
-----------------

Support for optional inequality constraints on individual components of the
solution vector has been added to CVODE and CVODES. See
:numref:`CVODE.Mathematics` and the description of in :numref:`CVODE.Usage.CC.optional_input` for
more details. Use of :c:func:`CVodeSetConstraints` requires the ``N_Vector``
operations :c:func:`N_VMinQuotient`, :c:func:`N_VConstMask`, and
:c:func:`N_VCompare` that were not previously required by CVODE and CVODES.

Fixed a problem with setting which would occur with some compilers (e.g.
armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA ``N_Vector`` library to from to better reflect that
we only support CUDA as a backend for RAJA currently.

Several changes were made to the build system:

   - CMake 3.1.3 is now the minimum required CMake version.
   - Deprecate the behavior of the CMake option and added the CMake option to select the integer size.
   - The native CMake FindMPI module is now used to locate an MPI installation.
   - If MPI is enabled and MPI compiler wrappers are not set, the build system will
     check if can compile MPI programs before trying to locate and use an MPI
     installation.
   - The previous options for setting MPI compiler wrappers and the executable for
     running MPI programs have been have been depreated. The new options that align
     with those used in native CMake FindMPI module are ``MPI_C_COMPILER``,
     ``MPO_CXX_COMPILER``, ``MPI_Fortran_COMPILER``, and ``MPIEXEC_EXECUTABLE``.
   - When a Fortran name-mangling scheme is needed (e.g., is ) the build system will
     infer the scheme from the Fortran compiler. If a Fortran compiler is not
     available or the inferred or default scheme needs to be overridden, the advanced
     options and can be used to manually set the name-mangling scheme and bypass
     trying to infer the scheme.
   - Parts of the main CMakeLists.txt file were moved to new files in the and
     directories to make the CMake configuration file structure more modular.

Changes in v3.1.2
-----------------

The changes in this minor release include the following:

-  Updated the minimum required version of CMake to 2.8.12 and enabled
   using rpath by default to locate shared libraries on OSX.
-  Fixed Windows specific problem where was not correctly defined when
   using 64-bit integers for the SUNDIALS index type. On Windows ``sunindextype`` is
   now defined as the MSVC basic type ``__int64``.
-  Added sparse SUNMatrix “Reallocate” routine to allow specification of
   the nonzero storage.
-  Updated the KLU SUNLinearSolver module to set constants for the two
   reinitialization types, and fixed a bug in the full reinitialization
   approach where the sparse SUNMatrix pointer would go out of scope on
   some architectures.
-  Updated the “ScaleAdd” and “ScaleAddI” implementations in the sparse
   SUNMatrix module to more optimally handle the case where the target
   matrix contained sufficient storage for the sum, but had the wrong
   sparsity pattern. The sum now occurs in-place, by performing the sum
   backwards in the existing storage. However, it is still more
   efficient if the user-supplied Jacobian routine allocates storage for
   the sum :math:`I+\gamma J` manually (with zero entries if needed).
-  Added the following examples from the usage notes page of the
   SUNDIALS website, and updated them to work with SUNDIALS 3.x:

   -  ``cvDisc_dns.c``, which demonstrates using CVODE with discontinuous solutions or RHS.
   -  ``cvRoberts_dns_negsol.c``, which illustrates the use of the RHS function return value to
      control unphysical negative concentrations.

-  Changed the LICENSE install path to `instdir/icnlude/sundials`.

Changes in v3.1.1
-----------------

The changes in this minor release include the following:

-  Fixed a minor bug in the cvSLdet routine, where a return was missing
   in the error check for three inconsistent roots.
-  Fixed a potential memory leak in the SPGMR and SPFGMR linear
   solvers: if “Initialize” was called multiple times then the solver
   memory was reallocated (without being freed).
-  Updated KLU ``SUNLinearSolver`` module to use a for the precision-specific
   solve function to be used (to avoid compiler warnings).
-  Added missing typecasts for some pointers (again, to avoid compiler
   warnings).
-  Bugfix in ``sunmatric_sparse.c`` where we had used instead of in one location.
-  Added missing ``#include <stio.h>`` in ``N_Vector`` and ``SUNMatrix`` header files.
-  Fixed an indexing bug in the CUDA ``N_Vector`` implementation of
   and revised the RAJA ``N_Vector`` implementation of :c:func:`N_VWrmsNormMask` to work with
   mask arrays using values other than zero or one. Replaced ``double`` with ``realtype`` in the
   RAJA vector test functions.
-  Fixed compilation issue with GCC 7.3.0 and Fortran programs that do
   not require a ``SUNMatrix`` or ``SUNLinearSolver`` module (e.g., iterative
   linear solvers or fixed-point iteration).

In addition to the changes above, minor corrections were also made to
the example programs, build system, and user documentation.

Changes in v3.1.0
-----------------

Added ``N_Vector`` print functions that write vector data to a specified
file (e.g., :c:func:`N_VPrintFile_Serial`).

Added ``make test`` and ``make test_install`` options to the build system for testing SUNDIALS after
building with and installing with respectively.

Changes in v3.0.0
-----------------

All interfaces to matrix structures and linear solvers have been reworked, and
all example programs have been updated. The goal of the redesign of these
interfaces was to provide more encapsulation and ease in interfacing custom
linear solvers and interoperability with linear solver libraries. Specific
changes include:

-  Added generic SUNMATRIX module with three provided implementations:
   dense, banded and sparse. These replicate previous SUNDIALS Dls and Sls
   matrix structures in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNMATRIX
   modules.

-  Added generic SUNLINEARSOLVER module with eleven provided
   implementations: dense, banded, LAPACK dense, LAPACK band, KLU, SuperLU_MT,
   SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG. These replicate previous SUNDIALS
   generic linear solvers in a single object-oriented API.

-  Added example problems demonstrating use of generic SUNLINEARSOLVER
   modules.

-  Expanded package-provided direct linear solver (Dls) interfaces and
   scaled, preconditioned, iterative linear solver (Spils) interfaces to utilize
   generic SUNMATRIX and SUNLINEARSOLVER objects.

-  Removed package-specific, linear solver-specific, solver modules
   (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since their functionality is
   entirely replicated by the generic Dls/Spils interfaces and
   SUNLINEARSOLVER/SUNMATRIX modules. The exception is CVDIAG, a diagonal
   approximate Jacobian solver available to CVODE and CVODES.

-  Converted all SUNDIALS example problems to utilize new generic
   SUNMATRIX and SUNLINEARSOLVER objects, along with updated Dls and Spils
   linear solver interfaces.

-  Added Spils interface routines to ARKode, CVODE, CVODES, IDA and IDAS
   to allow specification of a user-provided "JTSetup" routine. This change
   supports users who wish to set up data structures for the user-provided
   Jacobian-times-vector ("JTimes") routine, and where the cost of one JTSetup
   setup per Newton iteration can be amortized between multiple JTimes calls.

Two additional ``N_Vector`` implementations were added – one for CUDA and one
for RAJA vectors. These vectors are supplied to provide very basic support for
running on GPU architectures. Users are advised that these vectors both move all
data to the GPU device upon construction, and speedup will only be realized if
the user also conducts the right-hand-side function evaluation on the device. In
addition, these vectors assume the problem fits on one GPU. Further information
about RAJA, users are referred to th web site, https://software.llnl.gov/RAJA/.
These additions are accompanied by additions to various interface functions and
to user documentation.

All indices for data structures were updated to a new ``sunindextype`` that can
be configured to be a 32- or 64-bit integer data index type. ``sunindextype`` is
defined to be ``int32_t`` or ``int64_t`` when portable types are supported,
otherwise it is defined as ``int`` or ``long int``. The Fortran interfaces
continue to use for indices, except for their sparse matrix interface that now
uses the new . This new flexible capability for index types includes interfaces
to PETSc, hypre, SuperLU_MT, and KLU with either 32-bit or 64-bit capabilities
depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining ``booleantype``
values ``TRUE`` and ``FALSE`` have been changed to ``SUNTRUE`` and ``SUNFALSE``
respectively.

Temporary vectors were removed from preconditioner setup and solve routines for
all packages. It is assumed that all necessary data for user-provided
preconditioner operations will be allocated and stored in user-provided data
structures.

The file ``include/sundials_fconfig.h`` was added. This file contains SUNDIALS
type information for use in Fortran programs.

Added functions :c:func:`SUNDIALSGetVersion` and
:c:func:`SUNDIALSGetVersionNumber` to get SUNDIALS release version information
at runtime.

The build system was expanded to support many of the xSDK-compliant keys. The
xSDK is a movement in scientific software to provide a foundation for the rapid
and efficient production of high-quality, sustainable extreme-scale scientific
applications. More information can be found at, https://xsdk.info.

In addition, numerous changes were made to the build system. These include the
addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES`` CMake variables,
additional error checking during CMake configuration, minor bug fixes, and
renaming CMake options to enable/disable examples for greater clarity and an
added option to enable/disable Fortran 77 examples. These changes included
changing ``EXAMPLES_ENABLE`` to ``EXAMPLES_ENABLE_C``, changing ``CXX_ENABLE``
to ``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to ``EXAMPLES_ENABLE_F90``,
and adding an ``EXAMPLES_ENABLE_F77`` option.

A bug fix was made in :c:func:`CVodeFree` to call ``lfree`` unconditionally (if
non-NULL).

Corrections and additions were made to the examples, to installation-related
files, and to the user documentation.

Changes in v2.9.0
-----------------

Two additional ``N_Vector`` implementations were added – one for Hypre
(parallel) ParVector vectors, and one for PETSc vectors. These
additions are accompanied by additions to various interface functions
and to user documentation.

Each ``N_Vector`` module now includes a function, :c:func:`N_VGetVectorID`, that returns the
``N_Vector`` module name.

For each linear solver, the various solver performance counters are now
initialized to 0 in both the solver specification function and in solver
``linit`` function. This ensures that these solver counters are initialized upon
linear solver instantiation as well as at the beginning of the problem
solution.

In FCVODE, corrections were made to three Fortran interface
functions. Missing Fortran interface routines were added so that users
can supply the sparse Jacobian routine when using sparse direct solvers.

A memory leak was fixed in the banded preconditioner interface. In
addition, updates were done to return integers from linear solver and
preconditioner ’free’ functions.

The Krylov linear solver Bi-CGstab was enhanced by removing a redundant
dot product. Various additions and corrections were made to the
interfaces to the sparse solvers KLU and SuperLU_MT, including support
for CSR format when using KLU.

New examples were added for use of the OpenMP vector and for use of
sparse direct solvers from Fortran.

Minor corrections and additions were made to the CVODE solver, to
the Fortran interfaces, to the examples, to installation-related files,
and to the user documentation.

Changes in v2.8.0
-----------------

Two major additions were made to the linear system solvers that are available
for use with the CVODE solver. First, in the serial case, an interface to the
sparse direct solver KLU was added. Second, an interface to SuperLU_MT, the
multi-threaded version of SuperLU, was added as a thread-parallel sparse direct
solver option, to be used with the serial version of the ``N_Vector`` module. As
part of these additions, a sparse matrix (CSC format) structure was added to
CVODE.

Otherwise, only relatively minor modifications were made to the CVODE solver:

In ``cvRootFind``, a minor bug was corrected, where the input array was ignored,
and a line was added to break out of root-search loop if the initial interval
size is below the tolerance ``ttol``.

In ``CVLapackBand``, the line ``smu = MIN(N-1,mu+ml)`` was changed to to correct
an illegal input error for ``DGBTRF/DGBTRS``.

In order to eliminate or minimize the differences between the sources for
private functions in CVODE and CVODES, the names of 48 private functions were
changed from to , and a few other names were also changed.

Two minor bugs were fixed regarding the testing of input on the first call to –
one involving and one involving the initialization of ``*tret``.

In order to avoid possible name conflicts, the mathematical macro and function
names ``MIN``, ``MAX``, ``SQR``, ``RAbs``, ``RSqrt``, ``RExp``, ``RPowerI``, and
were changed to ``SUNMIN``, ``SUNMAX``, ``SUNSQR``, ``SUNRabs``, ``SUNRsqrt``,
``SUNRexp``, ``SUNRpowerI``, and ``SUNRPowerR`` respectively. These names occur
in both the solver and in various example programs.

The example program ``cvAdvDiff_diag_p`` was added to illustrate the use of in
parallel.

In the FCVODE optional input routines ``FCVSETIIN`` and ``FCVSETRIN``, the
optional fourth argument ``key_length`` was removed, with hardcoded key string
lengths passed to all tests.

In all FCVODE examples, integer declarations were revised so that those which
must match a C type ``long int`` are declared ``INTEGER*8``, and a comment was
added about the type match. All other integer declarations are just ``INTEGER``.
Corresponding minor corrections were made to the user guide.

Two new ``N_Vector`` modules have been added for thread-parallel computing
environments — one for OpenMP, denoted ``NVECTOR_OPENMP``, and one for Pthreads,
denoted ``NVECTOR_PTHREADS``.

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
function ``NewIntArray`` is replaced by a pair ``NewIntArray`` /
``NewLintArray``, for ``int`` and ``long int`` arrays, respectively.

A large number of minor errors have been fixed. Among these are the following:
In , the logic was changed to avoid a divide by zero. After the solver memory is
created, it is set to zero before being filled. In ``CVSetTqBDF`` each linear
solver interface function, the linear solver memory is freed on an error return,
and the function now includes a line setting to NULL the main memory pointer to
the linear solver memory. In the rootfinding functions ``CVRcheck1``/
``CVRcheck2``, when an exact zero is found, the array ``glo`` of :math:`g`
values at the left endpoint is adjusted, instead of shifting the :math:`t`
location slightly. In the installation files, we modified the treatment of the
macro SUNDIALS_USE_GENERIC_MATH, so that the parameter GENERIC_MATH_LIB is
either defined (with no value) or not defined.

Changes in v2.6.0
-----------------

Two new features were added in this release: (a) a new linear solver
module, based on BLAS and LAPACK for both dense and banded matrices, and
(b) an option to specify which direction of zero-crossing is to be
monitored while performing rootfinding.

The user interface has been further refined. Some of the API changes
involve: (a) a reorganization of all linear solver modules into two
families (besides the existing family of scaled preconditioned iterative
linear solvers, the direct solvers, including the new LAPACK-based ones,
were also organized into a *direct* family); (b) maintaining a single
pointer to user data, optionally specified through a -type function; and
(c) a general streamlining of the preconditioner modules distributed
with the solver.

Changes in v2.5.0
-----------------

The main changes in this release involve a rearrangement of the entire
:ref:`SUNDIALS source tree <CVODE.Organization>`. At the user interface
level, the main impact is in the mechanism of including SUNDIALS header files
which must now include the relative path (e.g. ``#include <cvode/cvode.h>``). Additional changes were made to
the build system: all exported header files are now installed in separate
subdirectories of the instaltion *include* directory.

The functions in the generic dense linear solver (``sundials_dense`` and
``sundials_smalldense``) were modified to work for rectangular :math:`m \times
n` matrices (:math:`m \le n`), while the factorization and solution functions
were renamed to ``DenseGETRF`` / ``denGETRF`` and ``DenseGETRS`` / ``denGETRS``,
respectively. The factorization and solution functions in the generic band
linear solver were renamed ``BandGBTRF`` and ``BandGBTRS``, respectively.

Changes in v2.4.0
-----------------

CVSPBCG and CVSPTFQMR modules have been added to interface with
the Scaled Preconditioned Bi-CGstab (SPBCG) and Scaled Preconditioned
Transpose-Free Quasi-Minimal Residual (SPTFQMR) linear solver modules,
respectively (for details see :numref:`CVODE.Usage.CC`). Corresponding additions were made
to the Fortran interface module FCVODE. At the same time, function type
names for Scaled Preconditioned Iterative Linear Solvers were added for
the user-supplied Jacobian-times-vector and preconditioner setup and
solve functions.

The deallocation functions now take as arguments the address of the
respective memory block pointer.

To reduce the possibility of conflicts, the names of all header files
have been changed by adding unique prefixes (``cvode_`` and ``sundials_``). When using the
default installation procedure, the header files are exported under
various subdirectories of the target directory. For more details see
:numref:`Installation`.

Changes in v2.3.0
-----------------

The user interface has been further refined. Several functions used for
setting optional inputs were combined into a single one. An optional
user-supplied routine for setting the error weight vector was added.
Additionally, to resolve potential variable scope issues, all SUNDIALS
solvers release user data right after its use. The build systems has
been further improved to make it more robust.

Changes in v2.2.1
-----------------

The changes in this minor SUNDIALS release affect only the build
system.

Changes in v2.2.0
-----------------

The major changes from the previous version involve a redesign of the
user interface across the entire SUNDIALS suite. We have eliminated
the mechanism of providing optional inputs and extracting optional
statistics from the solver through the and arrays. Instead, CVODE
now provides a set of routines (with prefix ``CVodeSet``) to change the default
values for various quantities controlling the solver and a set of
extraction routines (with prefix ``CVodeGet``) to extract statistics after return
from the main solver routine. Similarly, each linear solver module
provides its own set of - and -type routines. For more details see
:numref:`CVODE.Usage.CC.optional_input` and :numref:`CVODE.Usage.CC.optional_output`.

Additionally, the interfaces to several user-supplied routines (such as
those providing Jacobians and preconditioner information) were
simplified by reducing the number of arguments. The same information
that was previously accessible through such arguments can now be
obtained through -type functions.

The rootfinding feature was added, whereby the roots of a set of given
functions may be computed during the integration of the ODE system.

Installation of CVODE (and all of SUNDIALS) has been completely
redesigned and is now based on configure scripts.

.. _CVODE.Introduction.reading:

Reading this User Guide
=======================

This user guide is a combination of general usage instructions. Specific
example programs are provided as a separate document. We expect that
some readers will want to concentrate on the general instructions, while
others will refer mostly to the examples, and the organization is
intended to accommodate both styles.

There are different possible levels of usage of CVODE. The most
casual user, with a small IVP problem only, can get by with reading
:numref:`CVODE.Mathematics.ivp_sol`, then :numref:`CVODE.Usage.CC` through
:numref:`CVODE.Usage.CC.cvode` only, and looking at examples
in :cite:p:`cvode_ex`.

In a different direction, a more expert user with an IVP problem may
want to (a) use a package preconditioner
(:numref:`CVODE.Usage.CC.precond`), (b) supply his/her own Jacobian
or preconditioner routines
(:numref:`CVODE.Usage.CC.user_fct_sim.jacFn`), (c) do multiple runs of
problems of the same size (:numref:`CVODE.Usage.CC.reinit`), (d)
supply a new ``N_Vector`` module (:numref:`NVectors`), (e)
supply new ``SUNLinearSolver`` and/or ``SUNMatrix`` modules
(:numref:`SUNMatrix` and :numref:`SUNLinSol`),
or even (f) supply new ``SUNNonlinearSolver`` modules
(:numref:`SUNNonlinSol`).

The structure of this document is as follows:

-  In :numref:`CVODE.Mathematics`, we give short descriptions of the
   numerical methods implemented by CVODE for the solution of
   initial value problems for systems of ODEs, and continue with short
   descriptions of preconditioning
   (:numref:`CVODE.Mathematics.preconditioning`), stability limit
   detection (:numref:`CVODE.Mathematics.stablimit`), and rootfinding
   (:numref:`CVODE.Mathematics.rootfinding`).

-  The following chapter describes the structure of the SUNDIALS
   suite of solvers (:numref:`CVODE.Organization`) and the software
   organization of the CVODE solver (:numref:`CVODE.Organization.CVODE`).

-  :numref:`CVODE.Usage.CC` is the main usage document
   for CVODE for C applications. It includes a complete
   description of the user interface for the integration of ODE initial
   value problems.

-  In :numref:`SUNDIALS.Fortran`, we describe the use of
   CVODE with Fortran applications.

-  :numref:`NVectors` gives a brief overview of the
   generic ``N_Vector`` module shared among the various components of
   SUNDIALS, and details on the ``N_Vector`` implementations provided
   with SUNDIALS.

-  :numref:`SUNMatrix` gives a brief overview of
   the generic ``SUNMatrix`` module shared among the various components
   of SUNDIALS, and details on the ``SUNMatrix`` implementations
   provided with SUNDIALS: a dense implementation
   (:numref:`SUNMatrix.Dense`), a banded
   implementation (:numref:`SUNMatrix.Band`) and a
   sparse implementation
   (:numref:`SUNMatrix.Sparse`).

-  :numref:`SUNLinSol` gives a brief overview of
   the generic ``SUNLinearSolver`` module shared among the various components
   of SUNDIALS. This chapter contains details on the ``SUNLinearSolver``
   implementations provided with SUNDIALS. The chapter also contains
   details on the ``SUNLinearSolver`` implementations provided with
   SUNDIALS that interface with external linear solver libraries.

-  :numref:`SUNNonlinSol` describes the
   ``SUNNonlinearSolver`` API and nonlinear solver implementations shared
   among the various components of SUNDIALS.

-  Finally, in the appendices, we provide detailed instructions for the
   installation of CVODE, within the structure of SUNDIALS
   (:numref:`Installation`), as well as a list of all the
   constants used for input to and output from CVODE functions
   (:numref:`CVODE.Constants`).

Finally, the reader should be aware of the following notational
conventions in this user guide: program listings and identifiers (such
as :c:func:`CVodeInit`) within textual explanations are hyperlinked
to their definitions directly; fields
in C structures (such as *content*) appear in italics; and packages
or modules, such as CVLS, are written in all capitals.


SUNDIALS License and Notices
============================

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.


.. _CVODE.Introduction.Ack:

Acknowledgments
===============

We wish to acknowledge the contributions to previous versions of the
CVODE and PVODE codes and their user guides by Scott D.
Cohen :cite:p:`CoHi:94` and George D.
Byrne :cite:p:`ByHi:98`.
