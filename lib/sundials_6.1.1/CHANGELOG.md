# SUNDIALS Changelog

## Changes to SUNDIALS in release 6.1.1

Fixed exported `SUNDIALSConfig.cmake`.

Fixed Fortran interface to `MRIStepInnerStepper` and ``MRIStepCoupling`
structures and functions.

Added new Fortran example program,
`examples/arkode/F2003_serial/ark_kpr_mri_f2003.f90` demonstrating MRI
capabilities.

## Changes to SUNDIALS in release 6.1.0

Added new reduction implementations for the CUDA and HIP NVECTORs that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the `N_VSetKernelExecPolicy_Cuda` and
`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

`SUNDIALS::<lib>` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

``CMAKE_C_STANDARD`` is now set to 99 by default.

Fixed exported `SUNDIALSConfig.cmake` when profiling is enabled without Caliper.

Fixed `sundials_export.h` include in `sundials_config.h`.

Fixed memory leaks in the SUNLINSOL_SUPERLUMT linear solver.

## Changes to SUNDIALS in release 6.0.0

### SUNContext

SUNDIALS v6.0.0 introduces a new `SUNContext` object on which all other SUNDIALS
objects depend. As such, the constructors for all SUNDIALS packages, vectors,
matrices, linear solvers, nonlinear solvers, and memory helpers have been
updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call `SUNContext_Create` to create a context object with
before calling any other SUNDIALS library function, and then provide this object
to other SUNDIALS constructors. The context object has been introduced to allow
SUNDIALS to provide new features, such as the profiling/instrumentation also
introduced in this release, while maintaining thread-safety. See the
documentation section on the `SUNContext` for more details.

A script `upgrade-to-sundials-6-from-5.sh` has been provided with the release
(obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a `SUNCTX_PLACEHOLDER` argument to all of
the calls to SUNDIALS constructors that now require a `SUNContext` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

```
> ./upgrade-to-sundials-6-from-5.sh <files to update>
```

### SUNProfiler

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option `SUNDIALS_BUILD_WITH_PROFILING`. A built-in
profiler will be used by default, but the
[Caliper](https://github.com/LLNL/Caliper) library can also be used instead with
the CMake option `ENABLE_CALIPER`. See the documentation section on profiling
for more details.  **WARNING**: Profiling will impact performance, and should be
enabled judiciously.

### SUNMemoryHelper

The `SUNMemoryHelper` functions `Alloc`, `Dealloc`, and `Copy` have been updated
to accept an opaque handle as the last input. At a minimum, existing
`SUNMemoryHelper` implementations will need to update these functions to accept
the additional argument. Typically, this handle is the execution stream (e.g., a
CUDA/HIP stream or SYCL queue) for the operation. The CUDA, HIP, and SYCL
`SUNMemoryHelper` implementations have been updated accordingly. Additionally,
the constructor for the SYCL implementation has been updated to remove the SYCL
queue as an input.

### NVector

Two new optional vector operations, `N_VDotProdMultiLocal` and
`N_VDotProdMultiAllReduce`, have been added to support low-synchronization
methods for Anderson acceleration.

The CUDA, HIP, and SYCL execution policies have been moved from the `sundials`
namespace to the `sundials::cuda`, `sundials::hip`, and `sundials::sycl`
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The `Sundials` namespace used by the Trilinos Tpetra NVector has been replaced
with the `sundials::trilinos::nvector_tpetra` namespace.

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions `N_VCloneVectorArray_*` and `N_VDestroyVectorArray_*` have been
deprecated. The generic `N_VCloneVectorArray` and `N_VDestroyVectorArray`
functions should be used instead.

The previously deprecated constructor `N_VMakeWithManagedAllocator_Cuda` and
the function `N_VSetCudaStream_Cuda` have been removed and replaced with
`N_VNewWithMemHelp_Cuda` and `N_VSetKerrnelExecPolicy_Cuda` respectively.

The previously deprecated macros `PVEC_REAL_MPI_TYPE` and
`PVEC_INTEGER_MPI_TYPE` have been removed and replaced with
`MPI_SUNREALTYPE` and `MPI_SUNINDEXTYPE` respectively.

### SUNLinearSolver

The following previously deprecated functions have been removed

| Removed                   | Replaced with                    |
|:--------------------------|:---------------------------------|
| `SUNBandLinearSolver`     | `SUNLinSol_Band`                 |
| `SUNDenseLinearSolver`    | `SUNLinSol_Dense`                |
| `SUNKLU`                  | `SUNLinSol_KLU`                  |
| `SUNKLUReInit`            | `SUNLinSol_KLUReInit`            |
| `SUNKLUSetOrdering`       | `SUNLinSol_KLUSetOrdering`       |
| `SUNLapackBand`           | `SUNLinSol_LapackBand`           |
| `SUNLapackDense`          | `SUNLinSol_LapackDense`          |
| `SUNPCG`                  | `SUNLinSol_PCG`                  |
| `SUNPCGSetPrecType`       | `SUNLinSol_PCGSetPrecType`       |
| `SUNPCGSetMaxl`           | `SUNLinSol_PCGSetMaxl`           |
| `SUNSPBCGS`               | `SUNLinSol_SPBCGS`               |
| `SUNSPBCGSSetPrecType`    | `SUNLinSol_SPBCGSSetPrecType`    |
| `SUNSPBCGSSetMaxl`        | `SUNLinSol_SPBCGSSetMaxl`        |
| `SUNSPFGMR`               | `SUNLinSol_SPFGMR`               |
| `SUNSPFGMRSetPrecType`    | `SUNLinSol_SPFGMRSetPrecType`    |
| `SUNSPFGMRSetGSType`      | `SUNLinSol_SPFGMRSetGSType`      |
| `SUNSPFGMRSetMaxRestarts` | `SUNLinSol_SPFGMRSetMaxRestarts` |
| `SUNSPGMR`                | `SUNLinSol_SPGMR`                |
| `SUNSPGMRSetPrecType`     | `SUNLinSol_SPGMRSetPrecType`     |
| `SUNSPGMRSetGSType`       | `SUNLinSol_SPGMRSetGSType`       |
| `SUNSPGMRSetMaxRestarts`  | `SUNLinSol_SPGMRSetMaxRestarts`  |
| `SUNSPTFQMR`              | `SUNLinSol_SPTFQMR`              |
| `SUNSPTFQMRSetPrecType`   | `SUNLinSol_SPTFQMRSetPrecType`   |
| `SUNSPTFQMRSetMaxl`       | `SUNLinSol_SPTFQMRSetMaxl`       |
| `SUNSuperLUMT`            | `SUNLinSol_SuperLUMT`            |
| `SUNSuperLUMTSetOrdering` | `SUNLinSol_SuperLUMTSetOrdering` |

### Fortran Interfaces

The ARKODE, CVODE, IDA, and KINSOL Fortran 77 interfaces have been removed. See
the "SUNDIALS Fortran Interface" section in the user guides and the F2003
example programs for more details using the SUNDIALS Fortran 2003 module
interfaces.

### ARKODE

The ARKODE MRIStep module has been extended to support implicit-explicit (IMEX)
multirate infinitesimal generalized additive Runge-Kutta (MRI-GARK) methods. As
such, `MRIStepCreate` has been updated to include arguments for the slow
explicit and slow implicit ODE right-hand side functions. `MRIStepCreate` has
also been updated to require attaching an `MRIStepInnerStepper` for evolving the
fast time scale. `MRIStepReInit` has been similarly updated to take explicit
and implicit right-hand side functions as input. Codes using explicit or
implicit MRI methods will need to update `MRIStepCreate` and `MRIStepReInit`
calls to pass `NULL` for either the explicit or implicit right-hand side
function as appropriate. If ARKStep is used as the fast time scale integrator,
codes will need to call `ARKStepCreateMRIStepInnerStepper` to wrap the ARKStep
memory as an `MRIStepInnerStepper` object. Additionally, `MRIStepGetNumRhsEvals`
has been updated to return the number of slow implicit and explicit function
evaluations. The coupling table structure `MRIStepCouplingMem` and the
functions `MRIStepCoupling_Alloc` and `MRIStepCoupling_Create` have also
been updated to support IMEX-MRI-GARK methods.

The deprecated functions `MRIStepGetCurrentButcherTables` and
`MRIStepWriteButcher` and the utility functions `MRIStepSetTable` and
`MRIStepSetTableNum` have been removed. Users wishing to create an MRI-GARK
method from a Butcher table should use `MRIStepCoupling_MIStoMRI` to create
the corresponding MRI coupling table and attach it with `MRIStepSetCoupling`.

The implementation of solve-decoupled implicit MRI-GARK methods has been updated
to remove extraneous slow implicit function calls and reduce the memory
requirements.

Deprecated ARKODE nonlinear solver predictors: specification of the ARKStep
"bootstrap" or "minimum correction" predictors (options 4 and 5 from
`ARKStepSetPredictorMethod`), or MRIStep "bootstrap" predictor (option 4 from
`MRIStepSetPredictorMethod`), will output a deprecation warning message.
These options will be removed in a future release.

The previously deprecated functions `ARKStepSetMaxStepsBetweenLSet` and
`ARKStepSetMaxStepsBetweenJac` have been removed and replaced with
`ARKStepSetLSetupFrequency` and `ARKStepSetMaxStepsBetweenJac` respectively.

### CVODE

The previously deprecated function `CVodeSetMaxStepsBetweenJac` has been removed
and replaced with `CVodeSetJacEvalFrequency`.

### CVODES

Added a new function `CVodeGetLinSolveStats` to get the CVODES linear solver
statistics as a group.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODES after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODES statistics
throughout the simulation.

The previously deprecated function `CVodeSetMaxStepsBetweenJac` has been removed
and replaced with `CVodeSetJacEvalFrequency`.

### KINSOL

New orthogonalization methods were added for use within Anderson acceleration
in KINSOL. See the "Anderson Acceleration QR Factorization" subsection within
the mathematical considerations chapter of the user guide and the
`KINSetOrthAA` function documentation for more details.

### Deprecations

In addition to the deprecations noted elsewhere, many constants, types, and
functions have been renamed so that they are properly namespaced. The old names
have been deprecated and will be removed in SUNDIALS v7.0.0.

The following constants, macros, and  typedefs are now deprecated:

| Deprecated Name            | New Name                          |
|:---------------------------|:----------------------------------|
| `realtype`                 | `sunrealtype`                     |
| `booleantype`              | `sunbooleantype`                  |
| `RCONST`                   | `SUN_RCONST`                      |
| `BIG_REAL`                 | `SUN_BIG_REAL`                    |
| `SMALL_REAL`               | `SUN_SMALL_REAL`                  |
| `UNIT_ROUNDOFF`            | `SUN_UNIT_ROUNDOFF`               |
| `PREC_NONE`                | `SUN_PREC_NONE`                   |
| `PREC_LEFT`                | `SUN_PREC_LEFT`                   |
| `PREC_RIGHT`               | `SUN_PREC_RIGHT`                  |
| `PREC_BOTH`                | `SUN_PREC_BOTH`                   |
| `MODIFIED_GS`              | `SUN_MODIFIED_GS`                 |
| `CLASSICAL_GS`             | `SUN_CLASSICAL_GS`                |
| `ATimesFn`                 | `SUNATimesFn`                     |
| `PSetupFn`                 | `SUNPSetupFn`                     |
| `PSolveFn`                 | `SUNPSolveFn`                     |
| `DlsMat`                   | `SUNDlsMat`                       |
| `DENSE_COL`                | `SUNDLS_DENSE_COL`                |
| `DENSE_ELEM`               | `SUNDLS_DENSE_ELEM`               |
| `BAND_COL`                 | `SUNDLS_BAND_COL`                 |
| `BAND_COL_ELEM`            | `SUNDLS_BAND_COL_ELEM`            |
| `BAND_ELEM`                | `SUNDLS_BAND_ELEM`                |
| `SDIRK_2_1_2`              | `ARKODE_SDIRK_2_1_2`              |
| `BILLINGTON_3_3_2`         | `ARKODE_BILLINGTON_3_3_2`         |
| `TRBDF2_3_3_2`             | `ARKODE_TRBDF2_3_3_2`             |
| `KVAERNO_4_2_3`            | `ARKODE_KVAERNO_4_2_3`            |
| `ARK324L2SA_DIRK_4_2_3`    | `ARKODE_ARK324L2SA_DIRK_4_2_3`    |
| `CASH_5_2_4`               | `ARKODE_CASH_5_2_4`               |
| `CASH_5_3_4`               | `ARKODE_CASH_5_3_4`               |
| `SDIRK_5_3_4`              | `ARKODE_SDIRK_5_3_4`              |
| `KVAERNO_5_3_4`            | `ARKODE_KVAERNO_5_3_4`            |
| `ARK436L2SA_DIRK_6_3_4`    | `ARKODE_ARK436L2SA_DIRK_6_3_4`    |
| `KVAERNO_7_4_5`            | `ARKODE_KVAERNO_7_4_5`            |
| `ARK548L2SA_DIRK_8_4_5`    | `ARKODE_ARK548L2SA_DIRK_8_4_5`    |
| `ARK437L2SA_DIRK_7_3_4`    | `ARKODE_ARK437L2SA_DIRK_7_3_4`    |
| `ARK548L2SAb_DIRK_8_4_5`   | `ARKODE_ARK548L2SAb_DIRK_8_4_5`   |
| `MIN_DIRK_NUM`             | `ARKODE_MIN_DIRK_NUM`             |
| `MAX_DIRK_NUM`             | `ARKODE_MAX_DIRK_NUM`             |
| `MIS_KW3`                  | `ARKODE_MIS_KW3`                  |
| `MRI_GARK_ERK33a`          | `ARKODE_MRI_GARK_ERK33a`          |
| `MRI_GARK_ERK45a`          | `ARKODE_MRI_GARK_ERK45a`          |
| `MRI_GARK_IRK21a`          | `ARKODE_MRI_GARK_IRK21a`          |
| `MRI_GARK_ESDIRK34a`       | `ARKODE_MRI_GARK_ESDIRK34a`       |
| `MRI_GARK_ESDIRK46a`       | `ARKODE_MRI_GARK_ESDIRK46a`       |
| `IMEX_MRI_GARK3a`          | `ARKODE_IMEX_MRI_GARK3a`          |
| `IMEX_MRI_GARK3b`          | `ARKODE_IMEX_MRI_GARK3b`          |
| `IMEX_MRI_GARK4`           | `ARKODE_IMEX_MRI_GARK4`           |
| `MIN_MRI_NUM`              | `ARKODE_MIN_MRI_NUM`              |
| `MAX_MRI_NUM`              | `ARKODE_MAX_MRI_NUM`              |
| `DEFAULT_MRI_TABLE_3`      | `MRISTEP_DEFAULT_TABLE_3`         |
| `DEFAULT_EXPL_MRI_TABLE_3` | `MRISTEP_DEFAULT_EXPL_TABLE_3`    |
| `DEFAULT_EXPL_MRI_TABLE_4` | `MRISTEP_DEFAULT_EXPL_TABLE_4`    |
| `DEFAULT_IMPL_SD_TABLE_2`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_2` |
| `DEFAULT_IMPL_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_3` |
| `DEFAULT_IMPL_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_4` |
| `DEFAULT_IMEX_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_3` |
| `DEFAULT_IMEX_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_4` |
| `HEUN_EULER_2_1_2`         | `ARKODE_HEUN_EULER_2_1_2`         |
| `BOGACKI_SHAMPINE_4_2_3`   | `ARKODE_BOGACKI_SHAMPINE_4_2_3`   |
| `ARK324L2SA_ERK_4_2_3`     | `ARKODE_ARK324L2SA_ERK_4_2_3`     |
| `ZONNEVELD_5_3_4`          | `ARKODE_ZONNEVELD_5_3_4`          |
| `ARK436L2SA_ERK_6_3_4`     | `ARKODE_ARK436L2SA_ERK_6_3_4`     |
| `SAYFY_ABURUB_6_3_4`       | `ARKODE_SAYFY_ABURUB_6_3_4`       |
| `CASH_KARP_6_4_5`          | `ARKODE_CASH_KARP_6_4_5`          |
| `FEHLBERG_6_4_5`           | `ARKODE_FEHLBERG_6_4_5`           |
| `DORMAND_PRINCE_7_4_5`     | `ARKODE_DORMAND_PRINCE_7_4_5`     |
| `ARK548L2SA_ERK_8_4_5`     | `ARKODE_ARK548L2SA_ERK_8_4_5`     |
| `VERNER_8_5_6`             | `ARKODE_VERNER_8_5_6`             |
| `FEHLBERG_13_7_8`          | `ARKODE_FEHLBERG_13_7_8`          |
| `KNOTH_WOLKE_3_3`          | `ARKODE_KNOTH_WOLKE_3_3`          |
| `ARK437L2SA_ERK_7_3_4`     | `ARKODE_ARK437L2SA_ERK_7_3_4`     |
| `ARK548L2SAb_ERK_8_4_5`    | `ARKODE_ARK548L2SAb_ERK_8_4_5`    |
| `MIN_ERK_NUM`              | `ARKODE_MIN_ERK_NUM`              |
| `MAX_ERK_NUM`              | `ARKODE_MAX_ERK_NUM`              |
| `DEFAULT_ERK_2`            | `ARKSTEP_DEFAULT_ERK_2`           |
| `DEFAULT_ERK_3`            | `ARKSTEP_DEFAULT_ERK_3`           |
| `DEFAULT_ERK_4`            | `ARKSTEP_DEFAULT_ERK_4`           |
| `DEFAULT_ERK_5`            | `ARKSTEP_DEFAULT_ERK_5`           |
| `DEFAULT_ERK_6`            | `ARKSTEP_DEFAULT_ERK_6`           |
| `DEFAULT_ERK_8`            | `ARKSTEP_DEFAULT_ERK_8`           |
| `DEFAULT_DIRK_2`           | `ARKSTEP_DEFAULT_DIRK_2`          |
| `DEFAULT_DIRK_3`           | `ARKSTEP_DEFAULT_DIRK_3`          |
| `DEFAULT_DIRK_4`           | `ARKSTEP_DEFAULT_DIRK_4`          |
| `DEFAULT_DIRK_5`           | `ARKSTEP_DEFAULT_DIRK_5`          |
| `DEFAULT_ARK_ETABLE_3`     | `ARKSTEP_DEFAULT_ARK_ETABLE_3`    |
| `DEFAULT_ARK_ETABLE_4`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ETABLE_5`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ITABLE_3`     | `ARKSTEP_DEFAULT_ARK_ITABLE_3`    |
| `DEFAULT_ARK_ITABLE_4`     | `ARKSTEP_DEFAULT_ARK_ITABLE_4`    |
| `DEFAULT_ARK_ITABLE_5`     | `ARKSTEP_DEFAULT_ARK_ITABLE_5`    |
| `DEFAULT_ERK_2`            | `ERKSTEP_DEFAULT_2`               |
| `DEFAULT_ERK_3`            | `ERKSTEP_DEFAULT_3`               |
| `DEFAULT_ERK_4`            | `ERKSTEP_DEFAULT_4`               |
| `DEFAULT_ERK_5`            | `ERKSTEP_DEFAULT_5`               |
| `DEFAULT_ERK_6`            | `ERKSTEP_DEFAULT_6`               |
| `DEFAULT_ERK_8`            | `ERKSTEP_DEFAULT_8`               |

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

| Deprecated Name               | New Name                     |
|:------------------------------|:-----------------------------|
| `CVSpilsSetLinearSolver`      | `CVodeSetLinearSolver`       |
| `CVSpilsSetEpsLin`            | `CVodeSetEpsLin`             |
| `CVSpilsSetPreconditioner`    | `CVodeSetPreconditioner`     |
| `CVSpilsSetJacTimes`          | `CVodeSetJacTimes`           |
| `CVSpilsGetWorkSpace`         | `CVodeGetLinWorkSpace`       |
| `CVSpilsGetNumPrecEvals`      | `CVodeGetNumPrecEvals`       |
| `CVSpilsGetNumPrecSolves`     | `CVodeGetNumPrecSolves`      |
| `CVSpilsGetNumLinIters`       | `CVodeGetNumLinIters`        |
| `CVSpilsGetNumConvFails`      | `CVodeGetNumConvFails`       |
| `CVSpilsGetNumJTSetupEvals`   | `CVodeGetNumJTSetupEvals`    |
| `CVSpilsGetNumJtimesEvals`    | `CVodeGetNumJtimesEvals`     |
| `CVSpilsGetNumRhsEvals`       | `CVodeGetNumLinRhsEvals`     |
| `CVSpilsGetLastFlag`          | `CVodeGetLastLinFlag`        |
| `CVSpilsGetReturnFlagName`    | `CVodeGetLinReturnFlagName`  |
| `CVSpilsSetLinearSolverB`     | `CVodeSetLinearSolverB`      |
| `CVSpilsSetEpsLinB`           | `CVodeSetEpsLinB`            |
| `CVSpilsSetPreconditionerB`   | `CVodeSetPreconditionerB`    |
| `CVSpilsSetPreconditionerBS`  | `CVodeSetPreconditionerBS`   |
| `CVSpilsSetJacTimesB`         | `CVodeSetJacTimesB`          |
| `CVSpilsSetJacTimesBS`        | `CVodeSetJacTimesBS`         |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `CVDlsSetLinearSolverB`       | `CVodeSetLinearSolverB`      |
| `CVDlsSetJacFnB`              | `CVodeSetJacFnB`             |
| `CVDlsSetJacFnBS`             | `CVodeSetJacFnBS`            |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `KINDlsSetLinearSolver`       | `KINSetLinearSolver`         |
| `KINDlsSetJacFn`              | `KINSetJacFn`                |
| `KINDlsGetWorkSpace`          | `KINGetLinWorkSpace`         |
| `KINDlsGetNumJacEvals`        | `KINGetNumJacEvals`          |
| `KINDlsGetNumFuncEvals`       | `KINGetNumLinFuncEvals`      |
| `KINDlsGetLastFlag`           | `KINGetLastLinFlag`          |
| `KINDlsGetReturnFlagName`     | `KINGetLinReturnFlagName`    |
| `KINSpilsSetLinearSolver`     | `KINSetLinearSolver`         |
| `KINSpilsSetPreconditioner`   | `KINSetPreconditioner`       |
| `KINSpilsSetJacTimesVecFn`    | `KINSetJacTimesVecFn`        |
| `KINSpilsGetWorkSpace`        | `KINGetLinWorkSpace`         |
| `KINSpilsGetNumPrecEvals`     | `KINGetNumPrecEvals`         |
| `KINSpilsGetNumPrecSolves`    | `KINGetNumPrecSolves`        |
| `KINSpilsGetNumLinIters`      | `KINGetNumLinIters`          |
| `KINSpilsGetNumConvFails`     | `KINGetNumLinConvFails`      |
| `KINSpilsGetNumJtimesEvals`   | `KINGetNumJtimesEvals`       |
| `KINSpilsGetNumFuncEvals`     | `KINGetNumLinFuncEvals`      |
| `KINSpilsGetLastFlag`         | `KINGetLastLinFlag`          |
| `KINSpilsGetReturnFlagName`   | `KINGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolver`     | `IDASetLinearSolver`         |
| `IDASpilsSetPreconditioner`   | `IDASetPreconditioner`       |
| `IDASpilsSetJacTimes`         | `IDASetJacTimes`             |
| `IDASpilsSetEpsLin`           | `IDASetEpsLin`               |
| `IDASpilsSetIncrementFactor`  | `IDASetIncrementFactor`      |
| `IDASpilsGetWorkSpace`        | `IDAGetLinWorkSpace`         |
| `IDASpilsGetNumPrecEvals`     | `IDAGetNumPrecEvals`         |
| `IDASpilsGetNumPrecSolves`    | `IDAGetNumPrecSolves`        |
| `IDASpilsGetNumLinIters`      | `IDAGetNumLinIters`          |
| `IDASpilsGetNumConvFails`     | `IDAGetNumLinConvFails`      |
| `IDASpilsGetNumJTSetupEvals`  | `IDAGetNumJTSetupEvals`      |
| `IDASpilsGetNumJtimesEvals`   | `IDAGetNumJtimesEvals`       |
| `IDASpilsGetNumResEvals`      | `IDAGetNumLinResEvals`       |
| `IDASpilsGetLastFlag`         | `IDAGetLastLinFlag`          |
| `IDASpilsGetReturnFlagName`   | `IDAGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolverB`    | `IDASetLinearSolverB`        |
| `IDASpilsSetEpsLinB`          | `IDASetEpsLinB`              |
| `IDASpilsSetIncrementFactorB` | `IDASetIncrementFactorB`     |
| `IDASpilsSetPreconditionerB`  | `IDASetPreconditionerB`      |
| `IDASpilsSetPreconditionerBS` | `IDASetPreconditionerBS`     |
| `IDASpilsSetJacTimesB`        | `IDASetJacTimesB`            |
| `IDASpilsSetJacTimesBS`       | `IDASetJacTimesBS`           |
| `IDADlsSetLinearSolver`       | `IDASetLinearSolver`         |
| `IDADlsSetJacFn`              | `IDASetJacFn`                |
| `IDADlsGetWorkSpace`          | `IDAGetLinWorkSpace`         |
| `IDADlsGetNumJacEvals`        | `IDAGetNumJacEvals`          |
| `IDADlsGetNumResEvals`        | `IDAGetNumLinResEvals`       |
| `IDADlsGetLastFlag`           | `IDAGetLastLinFlag`          |
| `IDADlsGetReturnFlagName`     | `IDAGetLinReturnFlagName`    |
| `IDADlsSetLinearSolverB`      | `IDASetLinearSolverB`        |
| `IDADlsSetJacFnB`             | `IDASetJacFnB`               |
| `IDADlsSetJacFnBS`            | `IDASetJacFnBS`              |
| `DenseGETRF`                  | `SUNDlsMat_DenseGETRF`       |
| `DenseGETRS`                  | `SUNDlsMat_DenseGETRS`       |
| `denseGETRF`                  | `SUNDlsMat_denseGETRF`       |
| `denseGETRS`                  | `SUNDlsMat_denseGETRS`       |
| `DensePOTRF`                  | `SUNDlsMat_DensePOTRF`       |
| `DensePOTRS`                  | `SUNDlsMat_DensePOTRS`       |
| `densePOTRF`                  | `SUNDlsMat_densePOTRF`       |
| `densePOTRS`                  | `SUNDlsMat_densePOTRS`       |
| `DenseGEQRF`                  | `SUNDlsMat_DenseGEQRF`       |
| `DenseORMQR`                  | `SUNDlsMat_DenseORMQR`       |
| `denseGEQRF`                  | `SUNDlsMat_denseGEQRF`       |
| `denseORMQR`                  | `SUNDlsMat_denseORMQR`       |
| `DenseCopy`                   | `SUNDlsMat_DenseCopy`        |
| `denseCopy`                   | `SUNDlsMat_denseCopy`        |
| `DenseScale`                  | `SUNDlsMat_DenseScale`       |
| `denseScale`                  | `SUNDlsMat_denseScale`       |
| `denseAddIdentity`            | `SUNDlsMat_denseAddIdentity` |
| `DenseMatvec`                 | `SUNDlsMat_DenseMatvec`      |
| `denseMatvec`                 | `SUNDlsMat_denseMatvec`      |
| `BandGBTRF`                   | `SUNDlsMat_BandGBTRF`        |
| `bandGBTRF`                   | `SUNDlsMat_bandGBTRF`        |
| `BandGBTRS`                   | `SUNDlsMat_BandGBTRS`        |
| `bandGBTRS`                   | `SUNDlsMat_bandGBTRS`        |
| `BandCopy`                    | `SUNDlsMat_BandCopy`         |
| `bandCopy`                    | `SUNDlsMat_bandCopy`         |
| `BandScale`                   | `SUNDlsMat_BandScale`        |
| `bandScale`                   | `SUNDlsMat_bandScale`        |
| `bandAddIdentity`             | `SUNDlsMat_bandAddIdentity`  |
| `BandMatvec`                  | `SUNDlsMat_BandMatvec`       |
| `bandMatvec`                  | `SUNDlsMat_bandMatvec`       |
| `ModifiedGS`                  | `SUNModifiedGS`              |
| `ClassicalGS`                 | `SUNClassicalGS`             |
| `QRfact`                      | `SUNQRFact`                  |
| `QRsol`                       | `SUNQRsol`                   |
| `DlsMat_NewDenseMat`          | `SUNDlsMat_NewDenseMat`      |
| `DlsMat_NewBandMat`           | `SUNDlsMat_NewBandMat`       |
| `DestroyMat`                  | `SUNDlsMat_DestroyMat`       |
| `NewIntArray`                 | `SUNDlsMat_NewIntArray`      |
| `NewIndexArray`               | `SUNDlsMat_NewIndexArray`    |
| `NewRealArray`                | `SUNDlsMat_NewRealArray`     |
| `DestroyArray`                | `SUNDlsMat_DestroyArray`     |
| `AddIdentity`                 | `SUNDlsMat_AddIdentity`      |
| `SetToZero`                   | `SUNDlsMat_SetToZero`        |
| `PrintMat`                    | `SUNDlsMat_PrintMat`         |
| `newDenseMat`                 | `SUNDlsMat_newDenseMat`      |
| `newBandMat`                  | `SUNDlsMat_newBandMat`       |
| `destroyMat`                  | `SUNDlsMat_destroyMat`       |
| `newIntArray`                 | `SUNDlsMat_newIntArray`      |
| `newIndexArray`               | `SUNDlsMat_newIndexArray`    |
| `newRealArray`                | `SUNDlsMat_newRealArray`     |
| `destroyArray`                | `SUNDlsMat_destroyArray`     |

In addition, the entire `sundials_lapack.h` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

## Changes to SUNDIALS in release 5.8.0

The RAJA NVECTOR implementation has been updated to support the SYCL backend
in addition to the CUDA and HIP backend. Users can choose the backend when
configuring SUNDIALS by using the `SUNDIALS_RAJA_BACKENDS` CMake variable. This
module remains experimental and is subject to change from version to version.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix and the
linear solver support general dense linear systems as well as block diagonal
linear systems. This module is experimental and is subject to change from
version to version.

Added a new *optional* function to the SUNLinearSolver API,
`SUNLinSolSetZeroGuess`, to indicate that the next call to `SUNlinSolSolve` will
be made with a zero initial guess. SUNLinearSolver implementations that do not
use the `SUNLinSolNewEmpty` constructor will, at a minimum, need set the
`setzeroguess` function pointer in the linear solver `ops` structure to
`NULL`. The SUNDIALS iterative linear solver implementations have been updated
to leverage this new set function to remove one dot product per solve.

The time integrator packages (ARKODE, CVODE(S), and IDA(S)) all now support a
new "matrix-embedded" SUNLinearSolver type.  This type supports user-supplied
SUNLinearSolver implementations that set up and solve the specified linear
system at each linear solve call.  Any matrix-related data structures are held
internally to the linear solver itself, and are not provided by the SUNDIALS
package.

Added functions to ARKODE and CVODE(S) for supplying an alternative right-hand
side function and to IDA(S) for supplying an alternative residual for use within
nonlinear system function evaluations.

Support for user-defined inner (fast) integrators has been to the MRIStep module
in ARKODE. See the "MRIStep Custom Inner Steppers" section in the user guide for
more information on providing a user-defined integration method.

Added specialized fused HIP kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the `NVECTOR_HIP`
module. See the optional input function `CVodeSetUseIntegratorFusedKernels`
for more information. As with other SUNDIALS HIP features, this is
feature is experimental and may change from version to version.

New KINSOL options have been added to apply a constant damping factor in the
fixed point and Picard iterations (see `KINSetDamping`), to delay the start of
Anderson acceleration with the fixed point and Picard iterations (see
`KINSetDelayAA`), and to return the newest solution with the fixed point
iteration (see `KINSetReturnNewest`).

The installed SUNDIALSConfig.cmake file now supports the `COMPONENTS` option
to `find_package`. The exported targets no longer have IMPORTED_GLOBAL set.

A bug was fixed in `SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

A bug was fixed in the ARKODE stepper modules where the stop time may be passed
after resetting the integrator.

A bug was fixed in `IDASetJacTimesResFn` in IDAS where the supplied function was
used in the dense finite difference Jacobian computation rather than the finite
difference Jacobian-vector product approximation.

A bug was fixed in the KINSOL Picard iteration where the value of
`KINSetMaxSetupCalls` would be ignored.

## Changes to SUNDIALS in release 5.7.0

A new NVECTOR implementation based on the SYCL abstraction layer has been added
targeting Intel GPUs. At present the only SYCL compiler supported is the DPC++
(Intel oneAPI) compiler. See the SYCL NVECTOR section in the user guide for more
details. This module is considered experimental and is subject to major changes
even in minor releases.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the MAGMA linear algebra library. Both the matrix and the linear solver
support general dense linear systems as well as block diagonal linear systems,
and both are targeted at GPUs (AMD or NVIDIA).

## Changes to SUNDIALS in release 5.6.1

Fixed a bug in the SUNDIALS CMake which caused an error if the
`CMAKE_CXX_STANDARD` and `SUNDIALS_RAJA_BACKENDS` options were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

## Changes to SUNDIALS in release 5.6.0

A new NVECTOR implementation based on the AMD ROCm HIP platform has been added.
This vector can target NVIDIA or AMD GPUs. See HIP NVECTOR section in the user
guide for more details. This module is considered experimental and is subject to
change from version to version.

The RAJA NVECTOR implementation has been updated to support the HIP backend
in addition to the CUDA backend. Users can choose the backend when configuring
SUNDIALS by using the `SUNDIALS_RAJA_BACKENDS` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, `N_VGetDeviceArrayPointer`, was added to the N_Vector
API. This operation is useful for N_Vectors that utilize dual memory spaces,
e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations
no longer require the SUNDIALS CUDA N_Vector. Instead, they require that the
vector utilized provides the `N_VGetDeviceArrayPointer` operation, and that the
pointer returned by `N_VGetDeviceArrayPointer` is a valid CUDA device pointer.

## Changes to SUNDIALS in release 5.5.0

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes
should be fully backwards compatible for almost all users. SUNDIALS
now exports CMake targets and installs a `SUNDIALSConfig.cmake` file.

Added support for SuperLU DIST 6.3.0 or newer.

## Changes to SUNDIALS in release 5.4.0

Added full support for time-dependent mass matrices in ARKStep, and expanded
existing non-identity mass matrix infrastructure to support use of the
fixed point nonlinear solver.  Fixed bug for ERK method integration with
static mass matrices.

An interface between ARKStep and the XBraid multigrid reduction in time (MGRIT)
library has been added to enable parallel-in-time integration. See the ARKStep
documentation and examples for more details. This interface required the
addition of three new N_Vector operations to exchange vector data between
computational nodes, see `N_VBufSize`, `N_VBufPack`, and `N_VBufUnpack`. These
N_Vector operations are only used within the XBraid interface and need not be
implemented for any other context.

Updated the MRIStep time-stepping module in ARKODE to support higher-order
MRI-GARK methods [Sandu, SIAM J. Numer. Anal., 57, 2019], including methods that
involve solve-decoupled, diagonally-implicit treatment of the slow time scale.

A new API, `SUNMemoryHelper`, was added to support **GPU users** who have
complex memory management needs such as using memory pools. This is paired with
new constructors for the `NVECTOR_CUDA` and `NVECTOR_RAJA` modules that accept a
`SUNMemoryHelper` object. Refer to "The SUNMemoryHelper API", "NVECTOR CUDA" and
"NVECTOR RAJA" sections in the documentation for more information.

The `NVECTOR_RAJA` module has been updated to mirror the `NVECTOR_CUDA` module.
Notably, the update adds managed memory support to the `NVECTOR_RAJA` module.
Users of the module will need to update any calls to the `N_VMake_Raja` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

Added new `SetLSNormFactor` functions to CVODE(S), ARKODE, and IDA(S) to
to specify the factor for converting between integrator tolerances (WRMS norm)
and linear solver tolerances (L2 norm) i.e., `tol_L2 = nrmfac * tol_WRMS`.

Added new reset functions `ARKStepReset`, `ERKStepReset`, and
`MRIStepReset` to reset the stepper time and state vector to user-provided
values for continuing the integration from that point while retaining the
integration history. These function complement the reinitialization functions
`ARKStepReInit`, `ERKStepReInit`, and `MRIStepReInit` which reinitialize
the stepper so that the problem integration should resume as if started from
scratch.

Added new functions for advanced users providing a custom `SUNNonlinSolSysFn`.

The expected behavior of `SUNNonlinSolGetNumIters` and
`SUNNonlinSolGetNumConvFails` in the SUNNonlinearSolver API have been updated to
specify that they should return the number of nonlinear solver iterations and
convergence failures in the most recent solve respectively rather than the
cumulative number of iterations and failures across all solves respectively. The
API documentation and SUNDIALS provided SUNNonlinearSolver implementations and
have been updated accordingly. As before, the cumulative number of nonlinear
iterations and failures may be retrieved by calling the integrator provided get
functions.

**This change may cause a runtime error in existing user code**.
In IDAS and CVODES, the functions for forward integration with checkpointing
(`IDASolveF`, `CVodeF`) are now subject to a restriction on the number of time
steps allowed to reach the output time. This is the same restriction applied to
the `IDASolve` and `CVode` functions. The default maximum number of steps is
500, but this may be changed using the `<IDA|CVode>SetMaxNumSteps` function.
This change fixes a bug that could cause an infinite loop in the `IDASolveF`
and `CVodeF` and functions.

A minor inconsistency in CVODE(S) and a bug ARKODE when checking the Jacobian
evaluation frequency has been fixed. As a result codes using using a
non-default Jacobian update frequency through a call to
`CVodeSetMaxStepsBetweenJac` or `ARKStepSetMaxStepsBetweenJac` will need to
increase the provided value by 1 to achieve the same behavior as before. For
greater clarity the functions `CVodeSetMaxStepsBetweenJac`,
`ARKStepSetMaxStepsBetweenJac`, and `ARKStepSetMaxStepsBetweenLSet` have been
deprecated and replaced with `CVodeSetJacEvalFrequency`,
`ARKStepSetJacEvalFrequency`, and `ARKStepSetLSetupFrequency` respectively.
Additionally, the function `CVodeSetLSetupFrequency` has been added to CVODE(S)
to set the frequency of calls to the linear solver setup function.

The `NVECTOR_TRILINOS` module has been updated to work with Trilinos 12.18+.
This update changes the local ordinal type to always be an `int`.

Added support for CUDA v11.

## Changes to SUNDIALS in release 5.3.0

Fixed a bug in ARKODE where the prototypes for `ERKStepSetMinReduction` and
`ARKStepSetMinReduction` were not included in `arkode_erkstep.h` and
`arkode_arkstep.h` respectively.

Fixed a bug in ARKODE where inequality constraint checking would need to be
disabled and then re-enabled to update the inequality constraint values after
resizing a problem. Resizing a problem will now disable constraints and a call
to `ARKStepSetConstraints` or `ERKStepSetConstraints` is required to re-enable
constraint checking for the new problem size.

Fixed a bug in the iterative linear solver modules where an error is not
returned if the Atimes function is `NULL` or, if preconditioning is enabled, the
PSolve function is `NULL`.

Added specialized fused CUDA kernels to CVODE which may offer better
performance on smaller problems when using CVODE with the `NVECTOR_CUDA`
module. See the optional input function `CVodeSetUseIntegratorFusedKernels`
for more information. As with other SUNDIALS CUDA features, this is
feature is experimental and may change from version to version.

Added the ability to control the CUDA kernel launch parameters for the
`NVECTOR_CUDA` and `SUNMATRIX_CUSPARSE` modules. These modules remain
experimental and are subject to change from version to version.
In addition, the `NVECTOR_CUDA` kernels were rewritten to be more flexible.
Most users should see equivalent performance or some improvement, but a select
few may observe minor performance degradation with the default settings. Users
are encouraged to contact the SUNDIALS team about any performance changes
that they notice.

Added new capabilities for monitoring the solve phase in the
`SUNNONLINSOL_NEWTON` and `SUNNONLINSOL_FIXEDPOINT` modules, and the SUNDIALS
iterative linear solver modules. SUNDIALS must be built with the CMake option
`SUNDIALS_BUILD_WITH_MONITORING` to use these capabilities.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODE after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODE statistics
throughout the simulation.

Added a new function `CVodeGetLinSolveStats` to get the CVODE linear solver
statistics as a group.

Added optional set functions to provide an alternative ODE right-hand side
function (ARKODE and CVODE(S)), DAE residual function (IDA(S)), or nonlinear
system function (KINSOL) for use when computing Jacobian-vector products with
the internal difference quotient approximation.

Added support to CVODE for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user defined
projection function. This implementation is accompanied by additions to the
CVODE user documentation and examples.

## Changes to SUNDIALS in release 5.2.0

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set `CMAKE_Fortran_COMPILER` to `f2003`, `xlf2003`, or
`xlf2003_r`.

Fixed a bug in how ARKODE interfaces with a user-supplied, iterative, unscaled
linear solver. In this case, ARKODE adjusts the linear solver tolerance in an
attempt to account for the lack of support for left/right scaling matrices.
Previously, ARKODE computed this scaling factor using the error weight vector,
`ewt`; this fix changes that to the residual weight vector, `rwt`, that can
differ from `ewt` when solving problems with non-identity mass matrix.

Fixed a linkage bug affecting Windows users that stemmed from
dllimport/dllexport attribute missing on some SUNDIALS API functions.

Fixed a memory leak in CVODES and IDAS from not deallocating the `atolSmin0` and
`atolQSmin0` arrays.

Fixed a bug where a non-default value for the maximum allowed growth factor
after the first step would be ignored.

Functions were added to each of the time integration packages to enable or
disable the scaling applied to linear system solutions with matrix-based linear
solvers to account for lagged matrix information.

Added two new functions, `ARKStepSetMinReduction` and `ERKStepSetMinReduction`
to change the minimum allowed step size reduction factor after an error test
failure.

Added a new `SUNMatrix` implementation, `SUNMATRIX_CUSPARSE`, that interfaces to
the sparse matrix implementation from the NVIDIA cuSPARSE library. In addition,
the `SUNLINSOL_CUSOLVER_BATCHQR` linear solver has been updated to use this
matrix, therefore, users of this module will need to update their code. These
modules are still considered to be experimental, thus they are subject to
breaking changes even in minor releases.

Added a new "stiff" interpolation module to ARKODE, based on Lagrange polynomial
interpolation, that is accessible to each of the ARKStep, ERKStep and MRIStep
time-stepping modules. This module is designed to provide increased
interpolation accuracy when integrating stiff problems, as opposed to the ARKODE
standard Hermite interpolation module that can suffer when the IVP right-hand
side has large Lipschitz constant. While the Hermite module remains the default,
the new Lagrange module may be enabled using one of the routines
`ARKStepSetInterpolantType`, `ERKStepSetInterpolantType`, or
`MRIStepSetInterpolantType`. The serial example problem `ark_brusselator.c` has
been converted to use this Lagrange interpolation module. Created accompanying
routines `ARKStepSetInterpolantDegree`, `ARKStepSetInterpolantDegree` and
`ARKStepSetInterpolantDegree` to provide user control over these interpolating
polynomials. While the routines `ARKStepSetDenseOrder`, `ARKStepSetDenseOrder`
and `ARKStepSetDenseOrder` still exist, these have been deprecated and will be
removed in a future release.

## Changes to SUNDIALS in release 5.1.0

Added support for a user-supplied function to update the prediction for each
implicit stage solution in ARKStep.  If supplied, this routine will be called
*after* any existing ARKStep predictor algorithm completes, so that the
predictor may be modified by the user as desired.  The new user-supplied routine
has type `ARKStepStagePredictFn`, and may be set by calling
`ARKStepSetStagePredictFn`.

The MRIStep module has been updated to support attaching different user data
pointers to the inner and outer integrators. If applicable, user codes will
need to add a call to `ARKStepSetUserData` to attach their user data
pointer to the inner integrator memory as `MRIStepSetUserData` will
not set the pointer for both the inner and outer integrators. The MRIStep
examples have been updated to reflect this change.

Added support for damping when using Anderson acceleration in KINSOL. See the
mathematical considerations section of the user guide and the description of the
`KINSetDampingAA` function for more details.

Added support for damping to the `SUNNonlinearSolver_FixedPoint` module when
using Anderson acceleration. See the `SUNNonlinearSolver_FixedPoint` section in
the user guides and the description of the `SUNNonlinSolSetDamping_FixedPoint`
function for more details.

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables `PETSC_INCLUDES` and `PETSC_LIBRARIES` instead of `PETSC_DIR`.

Added a new build system option, `CUDA_ARCH`, to specify the CUDA architecture
to target.

Fixed a bug in the Fortran 2003 interfaces to the ARKODE Butcher table routines
and structure. This includes changing the `ARKodeButcherTable` type to be a
`type(c_ptr)` in Fortran.

Added two utility functions, `SUNDIALSFileOpen` and `SUNDIALSFileClose` for
creating/destroying file pointers. These are useful when using the Fortran 2003
interfaces.

## Changes to SUNDIALS in release 5.0.0

### Build System

Increased the minimum required CMake version to 3.5 for most SUNDIALS
configurations, and 3.10 when CUDA or OpenMP with device offloading are enabled.

The CMake option `BLAS_ENABLE` and the variable `BLAS_LIBRARIES` have been
removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
third party libraries that require linking to BLAS, the path to the BLAS
library should be included in the `_LIBRARIES` variable for the third party
library e.g., `SUPERLUDIST_LIBRARIES` when enabling SuperLU_DIST.

Fixed a bug in the build system that prevented the PThreads NVECTOR module from
being built.

### NVector

Two new functions were added to aid in creating custom NVECTOR objects. The
constructor `N_VNewEmpty` allocates an "empty" generic NVECTOR with the object's
content pointer and the function pointers in the operations structure
initialized to NULL. When used in the constructor for custom objects this
function will ease the introduction of any new optional operations to the
NVECTOR API by ensuring only required operations need to be set. Additionally,
the function `N_VCopyOps(w, v)` has been added to copy the operation function
pointers between vector objects. When used in clone routines for custom vector
objects these functions also will ease the introduction of any new optional
operations to the NVECTOR API by ensuring all operations are copied when cloning
objects.

Two new NVECTOR implementations, NVECTOR_MANYVECTOR and NVECTOR_MPIMANYVECTOR,
have been created to support flexible partitioning of solution data among
different processing elements (e.g., CPU + GPU) or for multi-physics problems
that couple distinct MPI-based simulations together (see the NVECTOR_MANYVECTOR
and NVECTOR_MPIMANYVECTOR sections in the user guides for more details). This
implementation is accompanied by additions to user documentation and SUNDIALS
examples.

An additional NVECTOR implementation, NVECTOR_MPIPLUSX, has been created to
support the MPI+X paradigm where X is a type of on-node parallelism (e.g.,
OpenMP, CUDA). The implementation is accompanied by additions to user
documentation and SUNDIALS examples.

One new required vector operation and ten new optional vector operations have
been added to the NVECTOR API. The new required operation, `N_VGetLength`,
returns the global length of an N_Vector. The optional operations have been
added to support the new NVECTOR_MPIMANYVECTOR implementation. The operation
`N_VGetCommunicator` must be implemented by subvectors that are combined to
create an NVECTOR_MPIMANYVECTOR, but is not used outside of this context. The
remaining nine operations are optional local reduction operations intended to
eliminate unnecessary latency when performing vector reduction operations
(norms, etc.) on distributed memory systems. The optional local reduction vector
operations are `N_VDotProdLocal`, `N_VMaxNormLocal`, `N_VMinLocal`,
`N_VL1NormLocal`, `N_VWSqrSumLocal`, `N_VWSqrSumMaskLocal`, `N_VInvTestLocal`,
`N_VConstrMaskLocal`, and `N_VMinQuotientLocal`. If an NVECTOR implementation
defines any of the local operations as NULL, then the NVECTOR_MPIMANYVECTOR will
call standard NVECTOR operations to complete the computation.

The `*_MPICuda` and `*_MPIRaja` functions have been removed from the
NVECTOR_CUDA and NVECTOR_RAJA implementations respectively. Accordingly, the
`nvector_mpicuda.h`, `nvector_mpiraja.h`, `libsundials_nvecmpicuda.lib`, and
`libsundials_nvecmpicudaraja.lib` files have been removed. Users should use the
NVECTOR_MPIPLUSX module in conjunction with the NVECTOR_CUDA or NVECTOR_RAJA
modules to replace the functionality. The necessary changes are minimal and
should require few code modifications.

Fixed a memory leak in the NVECTOR_PETSC clone function.

Made performance improvements to the CUDA NVECTOR. Users who utilize a
non-default stream should no longer see default stream synchronizations after
memory transfers.

Added a new constructor to the CUDA NVECTOR that allows a user to provide
custom allocate and free functions for the vector data array and internal
reduction buffer.

Added new Fortran 2003 interfaces for most NVECTOR modules. See NEVTOR section
in the user guides for more details on how to use the interfaces.

Added three new NVECTOR utility functions, `FN_VGetVecAtIndexVectorArray`,
`FN_VSetVecAtIndexVectorArray`, and `FN_VNewVectorArray`, for working with
`N_Vector` arrays when using the Fortran 2003 interfaces.

### SUNMatrix

Two new functions were added to aid in creating custom SUNMATRIX objects. The
constructor `SUNMatNewEmpty` allocates an "empty" generic SUNMATRIX with the
object's content pointer and the function pointers in the operations structure
initialized to NULL. When used in the constructor for custom objects this
function will ease the introduction of any new optional operations to the
SUNMATRIX API by ensuring only required operations need to be set. Additionally,
the function `SUNMatCopyOps(A, B)` has been added to copy the operation function
pointers between matrix objects. When used in clone routines for custom matrix
objects these functions also will ease the introduction of any new optional
operations to the SUNMATRIX API by ensuring all operations are copied when
cloning objects.

A new operation, `SUNMatMatvecSetup`, was added to the SUNMatrix API. Users
who have implemented custom SUNMatrix modules will need to at least update
their code to set the corresponding ops structure member, matvecsetup, to NULL.

The generic SUNMatrix API now defines error codes to be returned by SUNMatrix
operations. Operations which return an integer flag indiciating success/failure
may return different values than previously.

A new SUNMatrix (and SUNLinearSolver) implementation was added to facilitate
the use of the SuperLU_DIST library with SUNDIALS.

Added new Fortran 2003 interfaces for most SUNMATRIX modules. See SUNMATRIX
section in the user guides for more details on how to use the interfaces.

### SUNLinearSolver

A new function was added to aid in creating custom SUNLINEARSOLVER objects. The
constructor `SUNLinSolNewEmpty` allocates an "empty" generic SUNLINEARSOLVER
with the object's content pointer and the function pointers in the operations
structure initialized to NULL. When used in the constructor for custom objects
this function will ease the introduction of any new optional operations to the
SUNLINEARSOLVER API by ensuring only required operations need to be set.

The return type of the SUNLinearSolver API function `SUNLinSolLastFlag` has
changed from `long int` to `sunindextype` to be consistent with the type
used to store row indices in dense and banded linear solver modules.

Added a new optional operation to the SUNLINEARSOLVER API, `SUNLinSolGetID`,
that returns a `SUNLinearSolver_ID` for identifying the linear solver module.

The SUNLinearSolver API has been updated to make the initialize and setup
functions optional.

A new SUNLinearSolver (and SUNMatrix) implementation was added to facilitate
the use of the SuperLU_DIST library with SUNDIALS.

Added a new SUNLinearSolver implementation,
`SUNLinearSolver_cuSolverSp_batchQR`, which leverages the NVIDIA cuSOLVER sparse
batched QR method for efficiently solving block diagonal linear systems on
NVIDIA GPUs.

Added three new accessor functions to the SUNLinSol_KLU module,
`SUNLinSol_KLUGetSymbolic`, `SUNLinSol_KLUGetNumeric`, and
`SUNLinSol_KLUGetCommon`, to provide user access to the underlying
KLU solver structures.

Added new Fortran 2003 interfaces for most SUNLINEARSOLVER modules. See
SUNLINEARSOLVER section in the user guides for more details on how to use
the interfaces.

### SUNNonlinearSolver

A new function was added to aid in creating custom SUNNONLINEARSOLVER objects.
The constructor `SUNNonlinSolNewEmpty` allocates an "empty" generic
SUNNONLINEARSOLVER with the object's content pointer and the function pointers
in the operations structure initialized to NULL. When used in the constructor
for custom objects this function will ease the introduction of any new optional
operations to the SUNNONLINEARSOLVER API by ensuring only required operations
need to be set.

To facilitate the use of user supplied nonlinear solver convergence test
functions the `SUNNonlinSolSetConvTestFn` function in the SUNNonlinearSolver API
has been updated to take a `void*` data pointer as input. The supplied data
pointer will be passed to the nonlinear solver convergence test function on each
call.

The inputs values passed to the first two inputs of the `SUNNonlinSolSolve`
function in the SUNNONLINEARSOLVER have been changed to be the predicted
state and the initial guess for the correction to that state. Additionally,
the definitions of `SUNNonlinSolLSetupFn` and `SUNNonlinSolLSolveFn` in the
SUNNonlinearSolver API have been updated to remove unused input parameters.
For more information on the nonlinear system formulation and the API functions
see the SUNNONLINEARSOLVER chapter in the user guides.

Added a new `SUNNonlinearSolver` implementation, `SUNNonlinsol_PetscSNES`,
which interfaces to the PETSc SNES nonlinear solver API.

Added new Fortran 2003 interfaces for most SUNNONLINEARSOLVER modules. See
SUNNONLINEARSOLVER section in the user guides for more details on how to use
the interfaces.

### CVODE and CVODES

Fixed a bug in the CVODE and CVODES constraint handling where the step size
could be set below the minimum step size.

Fixed a bug in the CVODE and CVODES nonlinear solver interfaces where the norm
of the accumulated correction was not updated when using a non-default
convergence test function.

Fixed a bug in the CVODES `cvRescale` function where the loops to compute the
array of scalars for the fused vector scale operation stopped one iteration
early.

Fixed a bug in CVODES where CVodeF would return the wrong flag under certain
cirumstances.

Fixed a bug in CVODES where CVodeF would not return a root in NORMAL_STEP mode
if the root occurred after the desired output time.

Fixed a memeory leak in FCVODE when not using the default nonlinear solver.

Removed extraneous calls to `N_VMin` for simulations where the scalar valued
absolute tolerance, or all entries of the vector-valued absolute tolerance
array, are strictly positive. In this scenario CVODE and CVODES will remove
at least one global reduction per time step.

The CVLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type `SUNLINEARSOLVER_DIRECT`.

A new linear solver interface function, `CVLsLinSysFn`, was added as an
alternative method for evaluating the linear systems I - gamma J.

Added functions to get the current state and gamma value to CVODE and CVODES.
These functions may be useful to users who chose to provide their own nonlinear
solver implementation.

Added New Fortran 2003 interfaces to CVODE and CVODES were added. These new
interfaces were generated with SWIG-Fortran and provide a user an idiomatic
Fortran 2003 interface to most of the SUNDIALS C API. The existing CVODE F2003
interface, and all module implementations with existing Fortran 2003 interfaces
were updated accordingly. See the section "Using CVODE for Fortran
Applications" and "Using CVODES for Fortran Applications" in the appropriate
user guide for more details on how to use the interfaces.

### ARKODE

The MRIStep module has been updated to support explicit, implicit, or IMEX
methods as the fast integrator using the ARKStep module. As a result some
function signatures have been changed including MRIStepCreate which now
takes an ARKStep memory structure for the fast integration as an input.

Fixed a bug in the ARKStep time-stepping module in ARKODE that would result in
an infinite loop if the nonlinear solver failed to converge more than the
maximum allowed times during a single step.

Fixed a bug in ARKODE that would result in a "too much accuracy requested" error
when using fixed time step sizes with explicit methods in some cases.

Fixed a bug in ARKStep where the mass matrix linear solver setup function was
not called in the Matrix-free case.

Fixed a minor bug in ARKStep where an incorrect flag is reported when an
error occurs in the mass matrix setup or Jacobian-vector product setup
functions.

Fixed a memeory leak in FARKODE when not using the default nonlinear solver.

The reinitialization functions `ERKStepReInit`, `ARKStepReInit`, and
`MRIStepReInit` have been updated to retain the minimum and maxiumum step
size values from before reinitialization rather than resetting them to the
default values.

Removed extraneous calls to `N_VMin` for simulations where the scalar valued
absolute tolerance, or all entries of the vector-valued absolute tolerance
array, are strictly positive. In this scenario ARKODE steppers will remove
at least one global reduction per time step.

The ARKLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type `SUNLINEARSOLVER_DIRECT`.

A new linear solver interface function, `ARKLsLinSysFn`, was added as an
alternative method for evaluating the linear systems M - gamma J and
I - gamma J.

Added two new embedded ARK methods of orders 4 and 5 to ARKODE (from
Kennedy & Carpenter, Appl. Numer. Math., 136:183--205, 2019).

Support for optional inequality constraints on individual components of the
solution vector has been added the ARKODE ERKStep and ARKStep modules. See
the descriptions of `ERKStepSetConstraints` and `ARKStepSetConstraints` for
more details. Note that enabling constraint handling requires the NVECTOR
operations `N_VMinQuotient`, `N_VConstrMask`, and `N_VCompare` that were not
previously required by ARKODE.

Added functions to get the current state and gamma value to the ARKStep module.
These functions may be useful to users who chose to provide their own nonlinear
solver implementation.

Add two new 'Set' functions to MRIStep, `MRIStepSetPreInnerFn` and
`MRIStepSetPostInnerFn` for performing communication or memory
transfers needed before or after the inner integration.

Added new Fortran 2003 interfaces to all ARKODE stepper modules. These new
interfaces were generated with SWIG-Fortran and provide a user an idiomatic
Fortran 2003 interface to most of the SUNDIALS C API. See the section "Using
ARKODE for Fortran Applications" in the user guide for more details on how
to use the interfaces.

### IDA and IDAS

A bug was fixed in the IDA and IDAS linear solver interfaces where an incorrect
Jacobian-vector product increment was used with iterative solvers other than
SPGMR and SPFGMR.

Fixed a bug in IDAS where IDASolveF would return the wrong flag under certain
cirumstances.

Fixed a bug in IDAS where IDASolveF would not return a root in NORMAL_STEP mode
if the root occurred after the desired output time.

Fixed a bug the IDAS IDAQuadReInitB function where an incorrect memory structure
was passed to IDAQuadReInit.

Fixed a memeory leak in FIDA when not using the default nonlinear solver.

Removed extraneous calls to `N_VMin` for simulations where the scalar valued
absolute tolerance, or all entries of the vector-valued absolute tolerance
array, are strictly positive. In this scenario IDA and IDAS will remove
at least one global reduction per time step.

The IDALS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type SUNLINEARSOLVER_DIRECT.

Added new Fortran 2003 interfaces to IDA and IDAS. These new interfaces were
generated with SWIG-Fortran and provide a user an idiomatic Fortran 2003
interface to most of the SUNDIALS C API.  See the section "Using IDA for Fortran
Applications" and "Using IDAS for Fortran Applications" in the appropriate
user guide for more details on how to use the interfaces.

### KINSOL

Fixed a bug in the KINSOL linear solver interface where the auxiliary scalar
`sJpnorm` was not computed when necessary with the Picard iteration and the
auxiliary scalar `sFdotJp` was unnecessarily computed in some cases.

The KINLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type SUNLINEARSOLVER_DIRECT.

Added new Fortran 2003 interfaces to KINSOL. These new interfaces were
generated with SWIG-Fortran and provide a user an idiomatic Fortran 2003
interface to most of the SUNDIALS C API.  See the section "Using KINSOL for
Fortran Applications" for more details on how to use the interfaces.

## Changes to SUNDIALS in release 4.1.0

An additional N_Vector implementation was added for Tpetra vector from
Trilinos library to facilitate interoperability between SUNDIALS and Trilinos.
This implementation is accompanied by additions to user documentation and
SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The EXAMPLES_ENABLE_RAJA CMake option has been removed. The option
`EXAMPLES_ENABLE_CUDA` enables all examples that use CUDA including the RAJA
examples with a CUDA back end (if the RAJA NVECTOR is enabled).

The implementation header files (e.g. `arkode_impl.h`) are no longer installed.
This means users who are directly manipulating package memory structures will
need to update their code to use the package's public API.

Python is no longer required to run `make test` and `make test_install`.

Fixed a bug in `ARKodeButcherTable_Write` when printing a Butcher table
without an embedding.

## Changes to SUNDIALS in release 4.0.2

Added information on how to contribute to SUNDIALS and a contributing agreement.

Moved definitions of DLS and SPILS backwards compatibility functions to a source
file. The symbols are now included in the appropriate package library, e.g.
`libsundials_cvode.lib`.

## Changes to SUNDIALS in release 4.0.1

A bug in ARKODE where single precision builds would fail to compile has been
fixed.

## Changes to SUNDIALS in release 4.0.0

The direct and iterative linear solver interfaces in all SUNDIALS packages have
been merged into a single unified linear solver interface to support any valid
SUNLINSOL module. This includes the DIRECT and ITERATIVE types as well as the
new MATRIX_ITERATIVE type. Details regarding how SUNDIALS packages utilize
linear solvers of each type as well as discussion regarding intended use cases
for user-supplied SUNLINSOL implementations are included in the SUNLINSOL
chapter of the user guides. All example programs have been updated to use the
new unified interfaces.

The unified interface is very similar to the previous DLS and SPILS interfaces.
To minimize challenges in user migration to the unified linear solver interface,
the previous DLS and SPILS routines for all packages may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. Additionally, we note that Fortran users will need to enlarge
their iout array of optional integer outputs, and update the indices that they
query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided SUNLinSol
implementations have been updated to follow the naming convention SUNLinSol_*
where * is the name of the linear solver e.g., Dense, KLU, SPGMR, PCG, etc.
Solver-specific "set" routine names have been similarly standardized. To
minimize challenges in user migration to the new names, the previous routine
names may still be used; these will be deprecated in future releases, so we
recommend that users migrate to the new names soon. All example programs have
been updated to used the new naming convention.

The SUNBandMatrix constructor has been simplified to remove the storage upper
bandwidth argument.

SUNDIALS integrators (ARKODE, CVODE, CVODES, IDA, and IDAS) have been updated to
utilize generic nonlinear solver modules through the SUNNONLINSOL API. This API
will ease the addition of new nonlinear solver options and allow for external or
user-supplied nonlinear solvers. The SUNNONLINSOL API and provided SUNNONLINSOL
modules are described in a new user guide chapter and follow the same object
oriented design and implementation used by the NVECTOR, SUNMATRIX, and
SUNLINSOL modules. All integrator example programs have also been updated to
used the new nonlinear solver API.

Three fused vector operations and seven vector array operations have been added
to the NVECTOR API. These optional operations are disabled by default and may be
activated by calling vector specific routines after creating an NVECTOR. See the
NVECTOR chapter in the user guides for more information on the new operations.

Added a new NVECTOR (NVECTOR_OPENMPDEV) which leverages OpenMP 4.5+ device
offloading.

Multiple updates to the CUDA NVECTOR were made:

* Changed the `N_VMake_Cuda` function to take a host data pointer and a device
  data pointer instead of an `N_VectorContent_Cuda` object.

* Changed `N_VGetLength_Cuda` to return the global vector length instead of
  the local vector length.

* Added `N_VGetLocalLength_Cuda` to return the local vector length.

* Added `N_VGetMPIComm_Cuda` to return the MPI communicator used.

* Removed the accessor functions in the namespace suncudavec.

* Added the ability to set the `cudaStream_t` used for execution of the CUDA
  NVECTOR kernels. See the function `N_VSetCudaStreams_Cuda`.

* Added `N_VNewManaged_Cuda`, `N_VMakeManaged_Cuda`, and
  `N_VIsManagedMemory_Cuda` functions to accommodate using managed memory with
  the CUDA NVECTOR.

Multiple updates to the RAJA NVECTOR were made:

* Changed `N_VGetLength_Raja` to return the global vector length instead of
  the local vector length.

* Added `N_VGetLocalLength_Raja` to return the local vector length.

* Added `N_VGetMPIComm_Raja` to return the MPI communicator used.

* Removed the accessor functions in the namespace sunrajavec.

Two changes were made in the CVODE/CVODES/ARKODE initial step size algorithm:

  * Fixed an efficiency bug where an extra call to the RHS function was made.

  * Changed the behavior of the algorithm if the max-iterations case is hit.
    Before the algorithm would exit with the step size calculated on the
    penultimate iteration. Now it will exit with the step size calculated
    on the final iteration.

Fortran 2003 interfaces to CVODE, the fixed-point and Newton nonlinear solvers,
the dense, band, KLU, PCG, SPBCGS, SPFGMR, SPGMR, and SPTFQMR linear solvers,
and the serial, PThreads, and OpenMP NVECTORs have been added.

The ARKODE library has been entirely rewritten to support a modular approach to
one-step methods, which should allow for rapid research and development of novel
integration methods without affecting existing solver functionality.

A new ARKODE stepper, MRIStep, has been added for two rate explicit-explicit
multirate infinitesimal step methods.

ARKODE's dense output infrastructure has been improved to support higher-degree
Hermite polynomial interpolants (up to degree 5) over the last successful time
step.

## Changes to SUNDIALS in release 3.2.1

Fixed a bug in the CUDA NVECTOR where the `N_VInvTest` operation could write
beyond the allocated vector data.

Fixed library installation path for multiarch systems. This fix changes the
default library installation path to `CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR`
from `CMAKE_INSTALL_PREFIX/lib`. `CMAKE_INSTALL_LIBDIR` is automatically set,
but is available as a CMAKE option that can modified.

## Changes to SUNDIALS in release 3.2.0

Fixed problem with index types which would occur with some compilers (e.g.
armclang) that did not define `__STDC_VERSION__`. The fix includes a
depcrecation of the current behavior of the `SUNDIALS_INDEX_TYPE` CMake option.

Fixed a thread-safety issue in CVODES and IDAS when using adjoint sensitivity
analysis.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than one MPI
rank when using a GPU system. The vectors assume one GPU device per MPI rank.

Changed the name of the RAJA nvector library to `libsundials_nveccudaraja.lib`
from `libsundials_nvecraja.lib` to better reflect that we only support CUDA as a
backend for RAJA currently.

Increased CMake minimum version to 3.1.3

Add constraint handling feature to CVODE and CVODES.

Fixed a bug in IDAS where the saved residual value used in the nonlinear solve
for consistent initial conditions was passed as temporary workspace and could be
overwritten.

Several changes were made to the build system. If MPI is enabled and MPI
compiler wrappers are not set, the build system will check if
`CMAKE_<language>_COMPILER` can compile MPI programs before trying to locate and
use an MPI installation. The native CMake FindMPI module is now used to locate
an MPI installation. The options for setting MPI compiler wrappers and the
executable for running MPI programs have been updated to align with those in
native CMake FindMPI module. This included changing `MPI_MPICC` to
`MPI_C_COMPILER`, `MPI_MPICXX` to `MPI_CXX_COMPILER` combining `MPI_MPIF77` and
`MPI_MPIF90` to `MPI_Fortran_COMPILER`, and changing `MPI_RUN_COMMAND` to
`MPIEXEC_EXECUTABLE`. When a Fortran name-mangling scheme is needed (e.g.,
`LAPACK_ENABLE` is `ON`) the build system will infer the scheme from the Fortran
compiler. If a Fortran compiler is not available or the inferred or default
scheme needs to be overridden, the advanced options `SUNDIALS_F77_FUNC_CASE` and
`SUNDIALS_F77_FUNC_UNDERSCORES` can be used to manually set the name-mangling
scheme and bypass trying to infer the scheme. Additionally, parts of the main
`CMakeLists.txt` file were moved to new files in the src and example directories
to make the CMake configuration file structure more modular.

## Changes to SUNDIALS in release 3.1.2

Fixed Windows specific problem where `sunindextype` was not correctly defined
when using 64-bit integers. On Windows `sunindextype` is now defined as the MSVC
basic type `__int64`.

Changed LICENSE install path to `instdir/include/sundials`.

Updated the minimum required version of CMake to 2.8.12 and enabled using rpath
by default to locate shared libraries on OSX.

The misnamed function `CVSpilsSetJacTimesSetupFnBS` in cvodes has been
deprecated and replaced by `CVSpilsSetJacTimesBS`. The deprecated function
`CVSpilsSetJacTimesSetupFnBS` will be removed in the next major release.

Added and updated usage-notes examples from the SUNDIALS website to work with
SUNDIALS 3.x. The new examples are `cvode/cvDisc_dns.c`,
`cvode/cvRoberts_dns_negsol.c`, and `cvodes/cvsRoberts_FSA_dns_Switch.c`.

Added sparse SUNMatrix "Reallocate" routine to allow specification of the
nonzero storage.

Updated the KLU SUNLinearSolver module to set constants for the two
reinitialization types, and fixed a bug in the full reinitialization approach
where the sparse SUNMatrix pointer would go out of scope on some architectures.

Updated the "ScaleAdd" and "ScaleAddI" implementations in the sparse SUNMatrix
module to more optimally handle the case where the target matrix contained
sufficient storage for the sum, but had the wrong sparsity pattern. The sum now
occurs in-place, by performing the sum backwards in the existing storage.
However, it is still more efficient if the user-supplied Jacobian routine
allocates storage for the sum I + gamma J or M + gamma J manually (with zero
entries if needed).

## Changes to SUNDIALS in release 3.1.1

Fixed a minor bug in the CVODE and CVODES `cvSLdet` routine, where a return was
missing in the error check for three inconsistent roots.

Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers: if
"Initialize" was called multiple times then the solver memory was reallocated
(without being freed).

Fixed a minor bug in the `ARKReInit` routine, where a flag was incorrectly set
to indicate that the problem had been resized (instead of just re-initialized).

Fixed C++11 compiler errors/warnings about incompatible use of string literals.

Updated KLU SUNLinearSolver module to use a typedef for the precision-specific
solve function to be used (to avoid compiler warnings).

Added missing typecasts for some (`void*`) pointers to avoid compiler warnings.

Bugfix in `sunmatrix_sparse.c` where `int` was used instead of `sunindextype` in
one location.

Fixed a minor bug in `KINPrintInfo` where a case was missing for
`KIN_REPTD_SYSFUNC_ERR` leading to an undefined info message.

Added missing `#include <stdio.h>` in NVECTOR and SUNMATRIX header files.

Added missing prototypes for `ARKSpilsGetNumMTSetups` in ARKODE and
`IDASpilsGetNumJTSetupEvals` in IDA and IDAS.

Fixed an indexing bug in the CUDA NVECTOR implementation of `N_VWrmsNormMask`
and revised the RAJA NVECTOR implementation of `N_VWrmsNormMask` to work with
mask arrays using values other than zero or one. Replaced `double` with
`realtype` in the RAJA vector test functions.

Fixed compilation issue with GCC 7.3.0 and Fortran programs that do not require
a SUNMatrix or SUNLinearSolver module (e.g. iterative linear solvers, explicit
methods in ARKODE, functional iteration in CVODE, etc.).

## Changes to SUNDIALS in release 3.1.0

Added NVECTOR print functions that write vector data to a specified file (e.g.,
`N_VPrintFile_Serial`).

Added `make test` and `make test_install` options to the build system for
testing SUNDIALS after building with `make` and installing with `make install`
respectively.

Added "Changes in ..." (latest version) to all User Guides.

## Changes to SUNDIALS in release 3.0.0

Added new linear solver and matrix interfaces for all SUNDIALS packages and
updated the existing linear solver and matrix modules. The goal of the redesign
is to provide greater encapsulation and ease interfacing custom linear solvers
with linear solver libraries. Specific changes include:

 * Added generic SUNMATRIX module with three provided implementations:
   dense, banded and sparse.  These replicate previous SUNDIALS Dls and
   Sls matrix structures in a single object-oriented API.

 * Added example problems demonstrating use of generic SUNMATRIX modules.

 * Added generic SUNLINEARSOLVER module with eleven provided
   implementations: dense, banded, LAPACK dense, LAPACK band, KLU,
   SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG.  These replicate
   previous SUNDIALS generic linear solvers in a single object-oriented
   API.

 * Added example problems demonstrating use of generic SUNLINEARSOLVER
   modules.

 * Expanded package-provided direct linear solver (Dls) interfaces and
   scaled, preconditioned, iterative linear solver (Spils) interfaces
   to utilize generic SUNMATRIX and SUNLINEARSOLVER objects.

 * Removed package-specific, linear solver-specific, solver modules
   (e.g. CVDENSE, KINBAND, IDAKLU, ARKSPGMR) since their functionality
   is entirely replicated by the generic Dls/Spils interfaces and
   SUNLINEARSOLVER/SUNMATRIX modules.  The exception is CVDIAG, a
   diagonal approximate Jacobian solver available to CVODE and CVODES.

 * Converted all SUNDIALS example problems to utilize new generic
   SUNMATRIX and SUNLINEARSOLVER objects, along with updated Dls and
   Spils linear solver interfaces.

 * Added Spils interface routines to ARKODE, CVODE, CVODES, IDA and
   IDAS to allow specification of a user-provided "JTSetup" routine.
   This change supports users who wish to set up data structures for
   the user-provided Jacobian-times-vector ("JTimes") routine, and
   where the cost of one JTSetup setup per Newton iteration can be
   amortized between multiple JTimes calls.

Corresponding updates were made to all the example programs.

Two new NVECTOR modules added: for CUDA and RAJA support for GPU systems
(Information on RAJA: <https://software.llnl.gov/RAJA/> )
These vectors are supplied to provide very basic support for running
on GPU architectures.  Users are advised that these vectors both move all data
to the GPU device upon construction, and speedup will only be realized if the
user also conducts the right-hand-side function evaluation on the device.
In addition, these vectors assume the problem fits on one GPU.
For further information about RAJA, users are referred to the web site,
<https://software.llnl.gov/RAJA/.>

Addition of sunindextype option for 32-bit or 64-bit integer data index types
within all SUNDIALS structures

  * sunindextype is defined to be int32_t or int64_t when portable types are
    supported, otherwise it is defined as int or long int.

  * The Fortran interfaces continue to use `long int` for indices, except for
    their sparse matrix interface that now uses the new sunindextype.

  * Includes interfaces to PETSc, hypre, SuperLU_MT, and KLU with either 32-bit
    or 64-bit capabilities depending how the user configures SUNDIALS.

To avoid potential namespace conflicts, the macros defining booleantype
values TRUE and FALSE have been changed to SUNTRUE and SUNFALSE respectively.

Temporary vectors were removed from preconditioner setup and solve
routines for all packages.  It is assumed that all necessary data
for user-provided preconditioner operations will be allocated and
stored in user-provided data structures.

The file include/sundials\_fconfig.h was added.  This file contains
SUNDIALS type information for use in Fortran programs.

Added support for many xSDK-compliant build system keys
(Information on xSDK compliance: <https://xsdk.info/policies/> )
The xSDK is a movement in scientific software to provide a foundation for the
rapid and efficient production of high-quality,
sustainable extreme-scale scientific applications.  More information can
be found at <https://xsdk.info.>

Added functions SUNDIALSGetVersion and SUNDIALSGetVersionNumber to
get SUNDIALS release version information at runtime.

### Build System

Renamed CMake options to enable/disable examples for greater clarity
and added option to enable/disable Fortran 77 examples:

  * Changed `EXAMPLES_ENABLE` to `EXAMPLES_ENABLE_C`
  * Changed `CXX_ENABLE` to `EXAMPLES_ENABLE_CXX`
  * Changed `F90_ENABLE` to `EXAMPLES_ENABLE_F90`
  * Added `EXAMPLES_ENABLE_F77` option

Added separate `BLAS_ENABLE` and `BLAS_LIBRARIES` CMake variables

Fixed minor CMake bugs and included additional error checking during CMake
configuration

Corrections and additions to all User Guides.

Added "Changes in ..." (latest version) section to the introduction to in all
User Guides.

### ARKODE

Added comments to `arkode_butcher.c` regarding which methods should have
coefficients accurate enough for use in quad precision.

Fixed `RCONST` usage in `arkode_butcher.c`.

Fixed bug in `arkInitialSetup` to ensure the mass matrix vector product is
set up before the "msetup" routine is called.

Fixed ARKODE printf-related compiler warnings when building SUNDIALS
with extended precision.

### CVODE and CVODES

In `CVodeFree`, now call `lfree` unconditionally (if non-NULL).

### IDA and IDAS

Added missing prototype for `IDASetMaxBacksIC` in `ida.h` and `idas.h`.

### KINSOL

Corrected KINSOL fcmix name translation for `FKIN_SPFGMR`.

Renamed `KINLocalFn` and `KINCommFn` to `KINBBDLocalFn` and `KINBBDCommFn`
respectively in the BBD preconditioner module for consistency with other
SUNDIALS solvers.
