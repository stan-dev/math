.. ----------------------------------------------------------------
   Programmer(s): Daniel R. Reynolds @ SMU
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _ARKODE.Introduction:

Introduction
============

The ARKODE infrastructure provides adaptive-step time integration
modules for stiff, nonstiff and mixed stiff/nonstiff systems of
ordinary differential equations (ODEs).  ARKODE itself is structured
to support a wide range of one-step (but multi-stage) methods,
allowing for rapid development of parallel implementations of
state-of-the-art time integration methods.  At present, ARKODE is
packaged with two time-stepping modules, *ARKStep* and *ERKStep*.


*ARKStep* supports ODE systems posed in split, linearly-implicit form,

.. math::
   M(t)\, \dot{y} = f^E(t,y) + f^I(t,y),  \qquad y(t_0) = y_0,
   :label: ARKODE_ODE_split_linearly_implicit

where :math:`t` is the independent variable, :math:`y` is the set of
dependent variables (in :math:`\mathbb{R}^N`), :math:`M` is a
user-specified, nonsingular operator from :math:`\mathbb{R}^N` to
:math:`\mathbb{R}^N`, and the right-hand side function is partitioned
into up to two components:

- :math:`f^E(t,y)` contains the "nonstiff" time scale components to be
  integrated explicitly, and
- :math:`f^I(t,y)`  contains the "stiff" time scale components to be
  integrated implicitly.

Either of these operators may be disabled, allowing for fully
explicit, fully implicit, or combination implicit-explicit (ImEx) time
integration.

The algorithms used in ARKStep are adaptive- and fixed-step additive
Runge--Kutta methods. Such methods are defined through combining two
complementary Runge--Kutta methods: one explicit (ERK) and the other
diagonally implicit (DIRK).  Through appropriately partitioning the
ODE right-hand side into explicit and implicit components
:eq:`ARKODE_ODE_split_linearly_implicit`, such methods have the potential to
enable accurate and efficient time integration of stiff, nonstiff, and
mixed stiff/nonstiff systems of ordinary differential equations.  A
key feature allowing for high efficiency of these methods is that only
the components in :math:`f^I(t,y)` must be solved implicitly, allowing
for splittings tuned for use with optimal implicit solver algorithms.

This framework allows for significant freedom over the constitutive
methods used for each component, and ARKODE is packaged with a wide
array of built-in methods for use.  These built-in Butcher tables
include adaptive explicit methods of orders 2-8, adaptive implicit
methods of orders 2-5, and adaptive ImEx methods of orders 3-5.


*ERKStep* focuses specifically on problems posed in explicit form,

.. math::
   \dot{y} = f(t,y),  \qquad y(t_0) = y_0.
   :label: ARKODE_ODE_explicit

allowing for increased computational efficiency and memory savings.
The algorithms used in ERKStep are adaptive- and fixed-step explicit
Runge--Kutta methods.   As with ARKStep, the ERKStep module is packaged
with adaptive explicit methods of orders 2-8.


*MRIStep* focuses specifically on problems posed in additive form,

.. math::
   \dot{y} = f^E(t,y) + f^I(t,y) + f^F(t,y), \qquad y(t_0) = y_0.
   :label: ARKODE_ODE_two_rate

where here the right-hand side function is additively split into three
components:

* :math:`f^E(t,y)` contains the "slow-nonstiff" components of the system
  (this will be integrated using an explicit method and a large time step
  :math:`h^S`),

* :math:`f^I(t,y)` contains the "slow-stiff" components of the system
  (this will be integrated using an implicit method and a large time step
  :math:`h^S`), and

* :math:`f^F(t,y)` contains the "fast" components of the system (this will be
  integrated using a possibly different method than the slow time scale and a
  small time step :math:`h^F \ll h^S`).

For such problems, MRIStep provides fixed-step slow step multirate infinitesimal
step (MIS), multirate infinitesimal GARK (MRI-GARK), and implicit-explicit
MRI-GARK (IMEX-MRI-GARK) methods, allowing for evolution of the problem
:eq:`ARKODE_ODE_two_rate` using multirate methods having orders of accuracy 2-4.

For ARKStep or MRIStep problems that include nonzero implicit term
:math:`f^I(t,y)`, the resulting implicit system (assumed nonlinear, unless
specified otherwise) is solved approximately at each integration step, using a
SUNNonlinearSolver module, supplied either by the user or from the underlying
SUNDIALS infrastructure.  For nonlinear solver algorithms that internally
require a linear solver, ARKODE may use a variety of SUNLinearSolver modules
provided with SUNDIALS, or again may utilize a user-supplied module.




Changes from previous versions
--------------------------------

Changes in 5.0.0
^^^^^^^^^^^^^^^^

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

**ARKODE**

The MRIStep module has been extended to support implicit-explicit (ImEx)
multirate infinitesimal generalized additive Runge--Kutta (MRI-GARK) methods. As
such, :c:func:`MRIStepCreate` has been updated to include arguments for the slow
explicit and slow implicit ODE right-hand side functions.
:c:func:`MRIStepCreate` has also been updated to require attaching an
MRIStepInnerStepper for evolving the fast time scale. :c:func:`MRIStepReInit`
has been similarly updated to take explicit and implicit right-hand side
functions as input. Codes using explicit or implicit MRI methods will need to
update :c:func:`MRIStepCreate` and :c:func:`MRIStepReInit` calls to pass
``NULL`` for either the explicit or implicit right-hand side function as
appropriate. If ARKStep is used as the fast time scale integrator, codes will
need to call :c:func:`ARKStepCreateMRIStepInnerStepper` to wrap the ARKStep
memory as an MRIStepInnerStepper object. Additionally,
:c:func:`MRIStepGetNumRhsEvals` has been updated to return the number of slow
implicit and explicit function evaluations. The coupling table structure
:c:type:`MRIStepCouplingMem` and the functions :c:func:`MRIStepCoupling_Alloc`
and :c:func:`MRIStepCoupling_Create` have also been updated to support
IMEX-MRI-GARK methods.

The deprecated functions ``MRIStepGetCurrentButcherTables`` and
``MRIStepWriteButcher`` and the utility functions ``MRIStepSetTable`` and
``MRIStepSetTableNum`` have been removed. Users wishing to create an MRI-GARK
method from a Butcher table should use :c:func:`MRIStepCoupling_MIStoMRI` to
create the corresponding MRI coupling table and attach it with
:c:func:`MRIStepSetCoupling`.

The implementation of solve-decoupled implicit MRI-GARK methods has been updated
to remove extraneous slow implicit function calls and reduce the memory
requirements.

The previously deprecated functions ``ARKStepSetMaxStepsBetweenLSet`` and
``ARKStepSetMaxStepsBetweenJac`` have been removed and replaced with
:c:func:`ARKStepSetLSetupFrequency` and :c:func:`ARKStepSetMaxStepsBetweenJac`
respectively.

The ARKODE Fortran 77 interface has been removed. See :numref:`SUNDIALS.Fortran`
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
| ``SDIRK_2_1_2``              | ``ARKODE_SDIRK_2_1_2``              |
+------------------------------+-------------------------------------+
| ``BILLINGTON_3_3_2``         | ``ARKODE_BILLINGTON_3_3_2``         |
+------------------------------+-------------------------------------+
| ``TRBDF2_3_3_2``             | ``ARKODE_TRBDF2_3_3_2``             |
+------------------------------+-------------------------------------+
| ``KVAERNO_4_2_3``            | ``ARKODE_KVAERNO_4_2_3``            |
+------------------------------+-------------------------------------+
| ``ARK324L2SA_DIRK_4_2_3``    | ``ARKODE_ARK324L2SA_DIRK_4_2_3``    |
+------------------------------+-------------------------------------+
| ``CASH_5_2_4``               | ``ARKODE_CASH_5_2_4``               |
+------------------------------+-------------------------------------+
| ``CASH_5_3_4``               | ``ARKODE_CASH_5_3_4``               |
+------------------------------+-------------------------------------+
| ``SDIRK_5_3_4``              | ``ARKODE_SDIRK_5_3_4``              |
+------------------------------+-------------------------------------+
| ``KVAERNO_5_3_4``            | ``ARKODE_KVAERNO_5_3_4``            |
+------------------------------+-------------------------------------+
| ``ARK436L2SA_DIRK_6_3_4``    | ``ARKODE_ARK436L2SA_DIRK_6_3_4``    |
+------------------------------+-------------------------------------+
| ``KVAERNO_7_4_5``            | ``ARKODE_KVAERNO_7_4_5``            |
+------------------------------+-------------------------------------+
| ``ARK548L2SA_DIRK_8_4_5``    | ``ARKODE_ARK548L2SA_DIRK_8_4_5``    |
+------------------------------+-------------------------------------+
| ``ARK437L2SA_DIRK_7_3_4``    | ``ARKODE_ARK437L2SA_DIRK_7_3_4``    |
+------------------------------+-------------------------------------+
| ``ARK548L2SAb_DIRK_8_4_5``   | ``ARKODE_ARK548L2SAb_DIRK_8_4_5``   |
+------------------------------+-------------------------------------+
| ``MIN_DIRK_NUM``             | ``ARKODE_MIN_DIRK_NUM``             |
+------------------------------+-------------------------------------+
| ``MAX_DIRK_NUM``             | ``ARKODE_MAX_DIRK_NUM``             |
+------------------------------+-------------------------------------+
| ``MIS_KW3``                  | ``ARKODE_MIS_KW3``                  |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ERK33a``          | ``ARKODE_MRI_GARK_ERK33a``          |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ERK45a``          | ``ARKODE_MRI_GARK_ERK45a``          |
+------------------------------+-------------------------------------+
| ``MRI_GARK_IRK21a``          | ``ARKODE_MRI_GARK_IRK21a``          |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ESDIRK34a``       | ``ARKODE_MRI_GARK_ESDIRK34a``       |
+------------------------------+-------------------------------------+
| ``MRI_GARK_ESDIRK46a``       | ``ARKODE_MRI_GARK_ESDIRK46a``       |
+------------------------------+-------------------------------------+
| ``IMEX_MRI_GARK3a``          | ``ARKODE_IMEX_MRI_GARK3a``          |
+------------------------------+-------------------------------------+
| ``IMEX_MRI_GARK3b``          | ``ARKODE_IMEX_MRI_GARK3b``          |
+------------------------------+-------------------------------------+
| ``IMEX_MRI_GARK4``           | ``ARKODE_IMEX_MRI_GARK4``           |
+------------------------------+-------------------------------------+
| ``MIN_MRI_NUM``              | ``ARKODE_MIN_MRI_NUM``              |
+------------------------------+-------------------------------------+
| ``MAX_MRI_NUM``              | ``ARKODE_MAX_MRI_NUM``              |
+------------------------------+-------------------------------------+
| ``DEFAULT_MRI_TABLE_3``      | ``MRISTEP_DEFAULT_TABLE_3``         |
+------------------------------+-------------------------------------+
| ``DEFAULT_EXPL_MRI_TABLE_3`` | ``MRISTEP_DEFAULT_EXPL_TABLE_3``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_EXPL_MRI_TABLE_4`` | ``MRISTEP_DEFAULT_EXPL_TABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMPL_SD_TABLE_2``  | ``MRISTEP_DEFAULT_IMPL_SD_TABLE_2`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMPL_SD_TABLE_3``  | ``MRISTEP_DEFAULT_IMPL_SD_TABLE_3`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMPL_SD_TABLE_4``  | ``MRISTEP_DEFAULT_IMPL_SD_TABLE_4`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMEX_SD_TABLE_3``  | ``MRISTEP_DEFAULT_IMEX_SD_TABLE_3`` |
+------------------------------+-------------------------------------+
| ``DEFAULT_IMEX_SD_TABLE_4``  | ``MRISTEP_DEFAULT_IMEX_SD_TABLE_4`` |
+------------------------------+-------------------------------------+
| ``HEUN_EULER_2_1_2``         | ``ARKODE_HEUN_EULER_2_1_2``         |
+------------------------------+-------------------------------------+
| ``BOGACKI_SHAMPINE_4_2_3``   | ``ARKODE_BOGACKI_SHAMPINE_4_2_3``   |
+------------------------------+-------------------------------------+
| ``ARK324L2SA_ERK_4_2_3``     | ``ARKODE_ARK324L2SA_ERK_4_2_3``     |
+------------------------------+-------------------------------------+
| ``ZONNEVELD_5_3_4``          | ``ARKODE_ZONNEVELD_5_3_4``          |
+------------------------------+-------------------------------------+
| ``ARK436L2SA_ERK_6_3_4``     | ``ARKODE_ARK436L2SA_ERK_6_3_4``     |
+------------------------------+-------------------------------------+
| ``SAYFY_ABURUB_6_3_4``       | ``ARKODE_SAYFY_ABURUB_6_3_4``       |
+------------------------------+-------------------------------------+
| ``CASH_KARP_6_4_5``          | ``ARKODE_CASH_KARP_6_4_5``          |
+------------------------------+-------------------------------------+
| ``FEHLBERG_6_4_5``           | ``ARKODE_FEHLBERG_6_4_5``           |
+------------------------------+-------------------------------------+
| ``DORMAND_PRINCE_7_4_5``     | ``ARKODE_DORMAND_PRINCE_7_4_5``     |
+------------------------------+-------------------------------------+
| ``ARK548L2SA_ERK_8_4_5``     | ``ARKODE_ARK548L2SA_ERK_8_4_5``     |
+------------------------------+-------------------------------------+
| ``VERNER_8_5_6``             | ``ARKODE_VERNER_8_5_6``             |
+------------------------------+-------------------------------------+
| ``FEHLBERG_13_7_8``          | ``ARKODE_FEHLBERG_13_7_8``          |
+------------------------------+-------------------------------------+
| ``KNOTH_WOLKE_3_3``          | ``ARKODE_KNOTH_WOLKE_3_3``          |
+------------------------------+-------------------------------------+
| ``ARK437L2SA_ERK_7_3_4``     | ``ARKODE_ARK437L2SA_ERK_7_3_4``     |
+------------------------------+-------------------------------------+
| ``ARK548L2SAb_ERK_8_4_5``    | ``ARKODE_ARK548L2SAb_ERK_8_4_5``    |
+------------------------------+-------------------------------------+
| ``MIN_ERK_NUM``              | ``ARKODE_MIN_ERK_NUM``              |
+------------------------------+-------------------------------------+
| ``MAX_ERK_NUM``              | ``ARKODE_MAX_ERK_NUM``              |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_2``            | ``ARKSTEP_DEFAULT_ERK_2``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_3``            | ``ARKSTEP_DEFAULT_ERK_3``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_4``            | ``ARKSTEP_DEFAULT_ERK_4``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_5``            | ``ARKSTEP_DEFAULT_ERK_5``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_6``            | ``ARKSTEP_DEFAULT_ERK_6``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_8``            | ``ARKSTEP_DEFAULT_ERK_8``           |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_2``           | ``ARKSTEP_DEFAULT_DIRK_2``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_3``           | ``ARKSTEP_DEFAULT_DIRK_3``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_4``           | ``ARKSTEP_DEFAULT_DIRK_4``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_DIRK_5``           | ``ARKSTEP_DEFAULT_DIRK_5``          |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ETABLE_3``     | ``ARKSTEP_DEFAULT_ARK_ETABLE_3``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ETABLE_4``     | ``ARKSTEP_DEFAULT_ARK_ETABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ETABLE_5``     | ``ARKSTEP_DEFAULT_ARK_ETABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ITABLE_3``     | ``ARKSTEP_DEFAULT_ARK_ITABLE_3``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ITABLE_4``     | ``ARKSTEP_DEFAULT_ARK_ITABLE_4``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ARK_ITABLE_5``     | ``ARKSTEP_DEFAULT_ARK_ITABLE_5``    |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_2``            | ``ERKSTEP_DEFAULT_2``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_3``            | ``ERKSTEP_DEFAULT_3``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_4``            | ``ERKSTEP_DEFAULT_4``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_5``            | ``ERKSTEP_DEFAULT_5``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_6``            | ``ERKSTEP_DEFAULT_6``               |
+------------------------------+-------------------------------------+
| ``DEFAULT_ERK_8``            | ``ERKSTEP_DEFAULT_8``               |
+------------------------------+-------------------------------------+

In addition, the following functions are now deprecated (compile-time warnings
will be thrown if supported by the compiler):

+---------------------------------+--------------------------------+
| Deprecated Name                 | New Name                       |
+=================================+================================+
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

Changes in 4.8.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The RAJA NVECTOR implementation has been updated to support the SYCL backend
in addition to the CUDA and HIP backend. Users can choose the backend when
configuring SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable.
This module remains experimental and is subject to change from version to
version.

A new SUNMatrix and SUNLinearSolver implementation were added to interface with
the Intel oneAPI Math Kernel Library (oneMKL). Both the matrix and the linear
solver support general dense linear systems as well as block diagonal linear
systems. See :numref:`SUNLinSol.OneMklDense` for more details. This module is
experimental and is subject to change from version to version.

Added a new *optional* function to the SUNLinearSolver API,
:c:func:`SUNLinSolSetZeroGuess`, to indicate that the next call to
:c:func:`SUNLinSolSolve` will be made with a zero initial guess. SUNLinearSolver
implementations that do not use the :c:func:`SUNLinSolNewEmpty` constructor
will, at a minimum, need set the ``setzeroguess`` function pointer in the linear
solver ``ops`` structure to ``NULL``. The SUNDIALS iterative linear solver
implementations have been updated to leverage this new set function to remove
one dot product per solve.

ARKODE now supports a new "matrix-embedded" SUNLinearSolver type.  This type
supports user-supplied SUNLinearSolver implementations that set up and solve
the specified linear system at each linear solve call.  Any matrix-related data
structures are held internally to the linear solver itself, and are not
provided by the SUNDIALS package.

Support for user-defined inner (fast) integrators has been to the MRIStep
module. See :numref:`ARKODE.Usage.MRIStep.CustomInnerStepper` for more information on providing
a user-defined integration method.

Added the functions :c:func:`ARKStepSetNlsRhsFn()` and
:c:func:`MRIStepSetNlsRhsFn()` to supply an alternative implicit right-hand side
function for use within nonlinear system function evaluations.

The installed SUNDIALSConfig.cmake file now supports the ``COMPONENTS`` option
to ``find_package``. The exported targets no longer have IMPORTED_GLOBAL set.

A bug was fixed in :c:func:`SUNMatCopyOps` where the matrix-vector product setup
function pointer was not copied.

A bug was fixed in the SPBCGS and SPTFQMR solvers for the case where a non-zero
initial guess and a solution scaling vector are provided. This fix only impacts
codes using SPBCGS or SPTFQMR as standalone solvers as all SUNDIALS packages
utilize a zero initial guess.

A bug was fixed in the ARKODE stepper modules where the stop time may be passed
after resetting the integrator.


Changes in 4.7.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A new NVECTOR implementation based on the SYCL abstraction layer has been added
targeting Intel GPUs. At present the only SYCL compiler supported is the DPC++
(Intel oneAPI) compiler. See :numref:`NVectors.SYCL` for more details. This module
is considered experimental and is subject to major changes even in minor
releases.

A new SUNMatrix and SUNLinearSolver implementation were added to interface
with the MAGMA linear algebra library. Both the matrix and the linear solver
support general dense linear systems as well as block diagonal linear systems,
and both are targeted at GPUs (AMD or NVIDIA). See :numref:`SUNLinSol.MagmaDense`
for more details.

Changes in 4.6.1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed a bug in the SUNDIALS CMake which caused an error
if the CMAKE_CXX_STANDARD and SUNDIALS_RAJA_BACKENDS options
were not provided.

Fixed some compiler warnings when using the IBM XL compilers.

Changes in 4.6.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A new NVECTOR implementation based on the AMD ROCm HIP platform has been added.
This vector can target NVIDIA or AMD GPUs. See :numref:`NVectors.HIP` for more
details. This module is considered experimental and is subject to change from
version to version.

The RAJA NVECTOR implementation has been updated to support the HIP backend
in addition to the CUDA backend. Users can choose the backend when configuring
SUNDIALS by using the ``SUNDIALS_RAJA_BACKENDS`` CMake variable. This module
remains experimental and is subject to change from version to version.

A new optional operation, :c:func:`N_VGetDeviceArrayPointer`, was added to the
N_Vector API. This operation is useful for N_Vectors that utilize dual memory
spaces, e.g. the native SUNDIALS CUDA N_Vector.

The SUNMATRIX_CUSPARSE and SUNLINEARSOLVER_CUSOLVERSP_BATCHQR implementations
no longer require the SUNDIALS CUDA N_Vector. Instead, they require that the
vector utilized provides the :c:func:`N_VGetDeviceArrayPointer` operation, and
that the pointer returned by :c:func:`N_VGetDeviceArrayPointer` is a valid CUDA
device pointer.

Changes in v4.5.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Refactored the SUNDIALS build system. CMake 3.12.0 or newer is now required.
Users will likely see deprecation warnings, but otherwise the changes
should be fully backwards compatible for almost all users. SUNDIALS
now exports CMake targets and installs a SUNDIALSConfig.cmake file.

Added support for SuperLU DIST 6.3.0 or newer.

Changes in v4.4.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added full support for time-dependent mass matrices in ARKStep, and expanded
existing non-identity mass matrix infrastructure to support use of the
fixed point nonlinear solver. Fixed bug for ERK method integration with
static mass matrices.

An interface between ARKStep and the XBraid multigrid reduction in time (MGRIT)
library :cite:p:`xbraid` has been added to enable parallel-in-time integration. See the
:numref:`ARKODE.Usage.ARKStep.XBraid` section for more information and the example
codes in ``examples/arkode/CXX_xbraid``. This interface required the addition of
three new N_Vector operations to exchange vector data between computational
nodes, see :c:func:`N_VBufSize()`, :c:func:`N_VBufPack()`, and
:c:func:`N_VBufUnpack()`.  These N_Vector operations are only used within the
XBraid interface and need not be implemented for any other context.

Updated the MRIStep time-stepping module in ARKODE to support
higher-order MRI-GARK methods :cite:p:`Sandu:19`, including methods that
involve solve-decoupled, diagonally-implicit treatment of the
slow time scale.

Added the functions :c:func:`ARKStepSetLSNormFactor()`,
:c:func:`ARKStepSetMassLSNormFactor()`, and :c:func:`MRIStepSetLSNormFactor()`
to specify the factor for converting between integrator tolerances (WRMS norm)
and linear solver tolerances (L2 norm) i.e.,
``tol_L2 = nrmfac * tol_WRMS``.

Added new reset functions :c:func:`ARKStepReset()`, :c:func:`ERKStepReset()`,
and :c:func:`MRIStepReset()` to reset the stepper time and state vector to
user-provided values for continuing the integration from that point while
retaining the integration history. These function complement the
reinitialization functions :c:func:`ARKStepReInit()`, :c:func:`ERKStepReInit()`,
and :c:func:`MRIStepReInit()` which reinitialize the stepper so that the problem
integration should resume as if started from scratch.

Added new functions :c:func:`ARKStepComputeState`,
:c:func:`ARKStepGetNonlinearSystemData`, :c:func:`MRIStepComputeState`, and
:c:func:`MRIStepGetNonlinearSystemData` which advanced users might find useful
if providing a custom :c:func:`SUNNonlinSolSysFn`.

The expected behavior of :c:func:`SUNNonlinSolGetNumIters()` and
:c:func:`SUNNonlinSolGetNumConvFails()` in the SUNNonlinearSolver API have been
updated to specify that they should return the number of nonlinear solver
iterations and convergence failures in the most recent solve respectively rather
than the cumulative number of iterations and failures across all solves
respectively. The API documentation and SUNDIALS provided SUNNonlinearSolver
implementations have been updated accordingly. As before, the cumulative number
of nonlinear iterations may be retrieved by calling
:c:func:`ARKStepGetNumNonlinSolvIters()`, the cumulative number of failures with
:c:func:`ARKStepGetNumNonlinSolvConvFails()`, or both with
:c:func:`ARKStepGetNonlinSolvStats()`.

A minor bug in checking the Jacobian evaluation frequency has been fixed. As a
result codes using using a non-default Jacobian update frequency through a call
to :c:func:`ARKStepSetMaxStepsBetweenJac()` will need to increase the provided
value by 1 to achieve the same behavior as before. Additionally, for greater
clarity the functions :c:func:`ARKStepSetMaxStepsBetweenLSet()` and
:c:func:`ARKStepSetMaxStepsBetweenJac()` have been deprecated and replaced with
:c:func:`ARKStepSetLSetupFrequency()` and :c:func:`ARKStepSetJacEvalFrequency()`
respectively.

The ``NVECTOR_RAJA`` module has been updated to mirror the ``NVECTOR_CUDA`` module.
Notably, the update adds managed memory support to the ``NVECTOR_RAJA`` module.
Users of the module will need to update any calls to the ``N_VMake_Raja`` function
because that signature was changed. This module remains experimental and is
subject to change from version to version.

The ``NVECTOR_TRILINOS`` module has been updated to work with Trilinos 12.18+.
This update changes the local ordinal type to always be an ``int``.

Added support for CUDA v11.


Changes in v4.3.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed a bug in ARKODE where the prototypes for :c:func:`ERKStepSetMinReduction()`
and :c:func:`ARKStepSetMinReduction()` were not included in ``arkode_erkstep.h``
and ``arkode_arkstep.h`` respectively.

Fixed a bug where inequality constraint checking would need to be disabled and
then re-enabled to update the inequality constraint values after resizing a
problem. Resizing a problem will now disable constraints and a call to
:c:func:`ARKStepSetConstraints()` or :c:func:`ERKStepSetConstraints()` is
required to re-enable constraint checking for the new problem size.

Fixed a bug in the iterative linear solver modules where an error is not
returned if the Atimes function is ``NULL`` or, if preconditioning is enabled,
the PSolve function is ``NULL``.

Added the ability to control the CUDA kernel launch parameters for the
``NVECTOR_CUDA`` and ``SUNMATRIX_CUSPARSE`` modules. These modules remain
experimental and are subject to change from version to version.
In addition, the ``NVECTOR_CUDA`` kernels were rewritten to be more flexible.
Most users should see equivalent performance or some improvement, but a select
few may observe minor performance degradation with the default settings. Users
are encouraged to contact the SUNDIALS team about any perfomance changes
that they notice.

Added the optional function :c:func:`ARKStepSetJacTimesRhsFn()` to specify an
alternative implicit right-hand side function for computing Jacobian-vector
products with the internal difference quotient approximation.

Added new capabilities for monitoring the solve phase in the ``SUNNONLINSOL_NEWTON``
and ``SUNNONLINSOL_FIXEDPOINT`` modules, and the SUNDIALS iterative linear solver
modules. SUNDIALS must be built with the CMake option
``SUNDIALS_BUILD_WITH_MONITORING`` to use these capabilties.


Changes in v4.2.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed a build system bug related to the Fortran 2003 interfaces when using the
IBM XL compiler. When building the Fortran 2003 interfaces with an XL compiler
it is recommended to set ``CMAKE_Fortran_COMPILER`` to ``f2003``, ``xlf2003``,
or ``xlf2003_r``.

Fixed a bug in how ARKODE interfaces with a user-supplied, iterative, unscaled linear solver.
In this case, ARKODE adjusts the linear solver tolerance in an attempt to account for the
lack of support for left/right scaling matrices.  Previously, ARKODE computed this scaling
factor using the error weight vector, ``ewt``; this fix changes that to the residual weight vector,
``rwt``, that can differ from ``ewt`` when solving problems with non-identity mass matrix.

Fixed a similar bug in how ARKODE interfaces with scaled linear solvers when solving problems
with non-identity mass matrices.  Here, the left scaling matrix should correspond with ``rwt``
and the right scaling matrix with ``ewt``; these were reversed but are now correct.

Fixed a bug where a non-default value for the maximum allowed growth factor
after the first step would be ignored.

The function :c:func:`ARKStepSetLinearSolutionScaling()` was added to
enable or disable the scaling applied to linear system solutions with
matrix-based linear solvers to account for a lagged value of :math:`\gamma` in
the linear system matrix e.g., :math:`M - \gamma J` or :math:`I - \gamma J`.
Scaling is enabled by default when using a matrix-based linear solver.

Added two new functions, :c:func:`ARKStepSetMinReduction()` and
:c:func:`ERKStepSetMinReduction()`, to change the minimum allowed step size
reduction factor after an error test failure.

Added a new ``SUNMatrix`` implementation, :numref:`SUNMatrix.cuSparse`, that interfaces
to the sparse matrix implementation from the NVIDIA cuSPARSE library. In addition,
the :numref:`SUNLinSol.cuSolverSp` ``SUNLinearSolver`` has been updated to
use this matrix, as such, users of this module will need to update their code.
These modules are still considered to be experimental, thus they are subject to
breaking changes even in minor releases.

Added a new "stiff" interpolation module, based on Lagrange polynomial interpolation,
that is accessible to each of the ARKStep, ERKStep and MRIStep time-stepping modules.
This module is designed to provide increased interpolation accuracy when integrating
stiff problems, as opposed to the ARKODE-standard Hermite interpolation module that
can suffer when the IVP right-hand side has large Lipschitz constant.  While the
Hermite module remains the default, the new Lagrange module may be enabled using one
of the routines :c:func:`ARKStepSetInterpolantType()`, :c:func:`ERKStepSetInterpolantType()`,
or :c:func:`MRIStepSetInterpolantType()`.  The serial example problem ``ark_brusselator.c``
has been converted to use this Lagrange interpolation module.  Created accompanying routines
:c:func:`ARKStepSetInterpolantDegree()`, :c:func:`ARKStepSetInterpolantDegree()` and
:c:func:`ARKStepSetInterpolantDegree()` to provide user control over these
interpolating polynomials.  While the routines :c:func:`ARKStepSetDenseOrder()`,
:c:func:`ARKStepSetDenseOrder()` and :c:func:`ARKStepSetDenseOrder()` still exist,
these have been deprecated and will be removed in a future release.



Changes in v4.1.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed a build system bug related to finding LAPACK/BLAS.

Fixed a build system bug related to checking if the KLU library works.

Fixed a build system bug related to finding PETSc when using the CMake
variables ``PETSC_INCLUDES`` and ``PETSC_LIBRARIES`` instead of
``PETSC_DIR``.

Added a new build system option, ``CUDA_ARCH``, that can be used to specify
the CUDA architecture to compile for.

Fixed a bug in the Fortran 2003 interfaces to the ARKODE Butcher table routines and structure.
This includes changing the ``ARKodeButcherTable`` type to be a ``type(c_ptr)`` in Fortran.

Added two utility functions, ``SUNDIALSFileOpen`` and ``SUNDIALSFileClose``
for creating/destroying file pointers that are useful when using the Fortran
2003 interfaces.

Added support for a user-supplied function to update the prediction for each
implicit stage solution in ARKStep.  If supplied, this routine will be called
*after* any existing ARKStep predictor algorithm completes, so that the
predictor may be modified by the user as desired.  The new user-supplied routine
has type :c:type:`ARKStepStagePredictFn`, and may be set by calling
:c:func:`ARKStepSetStagePredictFn()`.

The MRIStep module has been updated to support attaching different user data
pointers to the inner and outer integrators. If applicable, user codes will
need to add a call to :c:func:`ARKStepSetUserData()` to attach their user data
pointer to the inner integrator memory as :c:func:`MRIStepSetUserData()` will
not set the pointer for both the inner and outer integrators. The MRIStep
examples have been updated to reflect this change.

Added support for constant damping to the ``SUNNonlinearSolver_FixedPoint``
module when using Anderson acceleration. See :numref:`SUNNonlinSol.FixedPoint.Math`
and the :c:func:`SUNNonlinSolSetDamping_FixedPoint()` for more details.

Changes in v4.0.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Build system changes**

Increased the minimum required CMake version to 3.5 for most SUNDIALS
configurations, and 3.10 when CUDA or OpenMP with device offloading are enabled.

The CMake option ``BLAS_ENABLE`` and the variable ``BLAS_LIBRARIES`` have been
removed to simplify builds as SUNDIALS packages do not use BLAS directly. For
third party libraries that require linking to BLAS, the path to the BLAS
library should be included in the ``_LIBRARIES`` variable for the third party
library e.g., ``SUPERLUDIST_LIBRARIES`` when enabling SuperLU_DIST.

Fixed a bug in the build system that prevented the PThreads NVECTOR module from
being built.

**NVECTOR module changes**

Two new functions were added to aid in creating custom NVECTOR objects. The
constructor :c:func:`N_VNewEmpty` allocates an "empty" generic NVECTOR with
the object's content pointer and the function pointers in the operations
structure initialized to ``NULL``. When used in the constructor for custom
objects this function will ease the introduction of any new optional operations
to the NVECTOR API by ensuring only required operations need to be set.
Additionally, the function :c:func:`N_VCopyOps()` has been added to copy the
operation function pointers between vector objects. When used in clone routines
for custom vector objects these functions also will ease the introduction of
any new optional operations to the NVECTOR API by ensuring all operations
are copied when cloning objects.

Two new NVECTOR implementations, NVECTOR_MANYVECTOR and
NVECTOR_MPIMANYVECTOR, have been created to support flexible partitioning
of solution data among different processing elements (e.g., CPU + GPU) or for
multi-physics problems that couple distinct MPI-based simulations together. This
implementation is accompanied by additions to user documentation and SUNDIALS
examples.

One new required vector operation and ten new optional vector operations have
been added to the NVECTOR API. The new required operation, :c:func:`N_VGetLength()`,
returns the global length of an ``N_Vector``. The optional operations have
been added to support the new NVECTOR_MPIMANYVECTOR implementation. The
operation :c:func:`N_VGetCommunicator()` must be implemented by subvectors that are
combined to create an NVECTOR_MPIMANYVECTOR, but is not used outside of
this context. The remaining nine operations are optional local reduction
operations intended to eliminate unnecessary latency when performing vector
reduction operations (norms, etc.) on distributed memory systems. The optional
local reduction vector operations are
:c:func:`N_VDotProdLocal`,
:c:func:`N_VMaxNormLocal`,
:c:func:`N_VMinLocal`,
:c:func:`N_VL1NormLocal`,
:c:func:`N_VWSqrSumLocal`,
:c:func:`N_VWSqrSumMaskLocal`,
:c:func:`N_VInvTestLocal`,
:c:func:`N_VConstrMaskLocal`, and
:c:func:`N_VMinQuotientLocal`.
If an NVECTOR implementation defines any of the local operations as
``NULL``, then the NVECTOR_MPIMANYVECTOR will call standard NVECTOR
operations to complete the computation.

An additional NVECTOR implementation, NVECTOR_MPIPLUSX, has been created to
support the MPI+X paradigm where X is a type of on-node parallelism
(*e.g.*, OpenMP, CUDA). The implementation is accompanied by additions to
user documentation and SUNDIALS examples.

The ``*_MPICuda`` and ``*_MPIRaja`` functions have been removed from the
NVECTOR_CUDA and NVECTOR_RAJA implementations respectively. Accordingly, the
``nvector_mpicuda.h``, ``nvector_mpiraja.h``, ``libsundials_nvecmpicuda.lib``,
and ``libsundials_nvecmpicudaraja.lib`` files have been removed. Users should
use the NVECTOR_MPIPLUSX module coupled in conjunction with the NVECTOR_CUDA
or NVECTOR_RAJA modules to replace the functionality. The necessary changes are
minimal and should require few code modifications. See the programs in
``examples/ida/mpicuda`` and ``examples/ida/mpiraja`` for examples of how to
use the NVECTOR_MPIPLUSX module with the NVECTOR_CUDA and NVECTOR_RAJA modules
respectively.

Fixed a memory leak in the NVECTOR_PETSC module clone function.

Made performance improvements to the NVECTOR_CUDA module. Users who utilize a
non-default stream should no longer see default stream synchronizations
after memory transfers.

Added a new constructor to the NVECTOR_CUDA module that allows a user to provide
custom allocate and free functions for the vector data array and internal
reduction buffer.

Added new Fortran 2003 interfaces for most NVECTOR modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

Added three new NVECTOR utility functions,
:c:func:`N_VGetVecAtIndexVectorArray()`
:c:func:`N_VSetVecAtIndexVectorArray()`, and
:c:func:`N_VNewVectorArray`,
for working with ``N_Vector`` arrays when using the Fortran 2003 interfaces.

**SUNMatrix module changes**

Two new functions were added to aid in creating custom SUNMATRIX objects. The
constructor :c:func:`SUNMatNewEmpty` allocates an "empty" generic SUNMATRIX with
the object's content pointer and the function pointers in the operations
structure initialized to ``NULL``. When used in the constructor for custom
objects this function will ease the introduction of any new optional operations
to the SUNMATRIX API by ensuring only required operations need to be set.
Additionally, the function :c:func:`SUNMatCopyOps()` has been added to copy the
operation function pointers between matrix objects. When used in clone routines
for custom matrix objects these functions also will ease the introduction of any
new optional operations to the SUNMATRIX API by ensuring all operations are
copied when cloning objects.

A new operation, :c:func:`SUNMatMatvecSetup()`, was added to the SUNMATRIX API.
Users who have implemented custom SUNMATRIX modules will need to at least
update their code to set the corresponding ``ops`` structure member,
``matvecsetup``, to ``NULL``.

A new operation, :c:func:`SUNMatMatvecSetup()`, was added to the SUNMATRIX API
to perform any setup necessary for computing a matrix-vector product. This
operation is useful for SUNMATRIX implementations which need to prepare the
matrix itself, or communication structures before performing the matrix-vector
product. Users who have implemented custom SUNMATRIX modules will need to at
least update their code to set the corresponding ``ops`` structure member,
``matvecsetup``, to ``NULL``.

The generic SUNMATRIX API now defines error codes to be returned by
SUNMATRIX operations. Operations which return an integer flag indiciating
success/failure may return different values than previously.

A new SUNMATRIX (and SUNLINEARSOLVER) implementation was added to
facilitate the use of the SuperLU_DIST library with SUNDIALS.

Added new Fortran 2003 interfaces for most SUNMATRIX modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

**SUNLinearSolver module changes**

A new function was added to aid in creating custom SUNLINEARSOLVER objects.
The constructor :c:func:`SUNLinSolNewEmpty` allocates an "empty" generic
SUNLINEARSOLVER with the object's content pointer and the function pointers
in the operations structure initialized to ``NULL``. When used in the
constructor for custom objects this function will ease the introduction of any
new optional operations to the SUNLINEARSOLVER API by ensuring only required
operations need to be set.

The return type of the SUNLINEARSOLVER API function :c:func:`SUNLinSolLastFlag()`
has changed from ``long int`` to ``sunindextype`` to be consistent with the
type used to store row indices in dense and banded linear solver modules.

Added a new optional operation to the SUNLINEARSOLVER API,
:c:func:`SUNLinSolGetID`, that returns a ``SUNLinearSolver_ID`` for identifying
the linear solver module.

The SUNLINEARSOLVER API has been updated to make the initialize and setup
functions optional.

A new SUNLINEARSOLVER (and SUNMATRIX) implementation was added to
facilitate the use of the SuperLU_DIST library with SUNDIALS.

Added a new SUNLinearSolver implementation, ``SUNLinearSolver_cuSolverSp_batchQR``,
which leverages the NVIDIA cuSOLVER sparse batched QR method for efficiently
solving block diagonal linear systems on NVIDIA GPUs.


Added three new accessor functions to the SUNLinSol_KLU module,
:c:func:`SUNLinSol_KLUGetSymbolic()`, :c:func:`SUNLinSol_KLUGetNumeric()`, and
:c:func:`SUNLinSol_KLUGetCommon()`, to provide user access to the underlying
KLU solver structures.

Added new Fortran 2003 interfaces for most SUNLINEARSOLVER modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

**SUNNonlinearSolver module changes**

A new function was added to aid in creating custom SUNNONLINEARSOLVER
objects. The constructor :c:func:`SUNNonlinSolNewEmpty` allocates an "empty"
generic SUNNONLINEARSOLVER with the object's content pointer and the function
pointers in the operations structure initialized to ``NULL``. When used in the
constructor for custom objects this function will ease the introduction of any
new optional operations to the SUNNONLINEARSOLVER API by ensuring only
required operations need to be set.

To facilitate the use of user supplied nonlinear solver convergence test
functions the :c:func:`SUNNonlinSolSetConvTestFn()` function in the
SUNNONLINEARSOLVER API has been updated to take a ``void*`` data pointer as
input. The supplied data pointer will be passed to the nonlinear solver
convergence test function on each call.

The inputs values passed to the first two inputs of the :c:func:`SUNNonlinSolSolve()`
function in the SUNNONLINEARSOLVER have been changed to be the predicted
state and the initial guess for the correction to that state. Additionally,
the definitions of :c:type:`SUNNonlinSolLSetupFn` and :c:type:`SUNNonlinSolLSolveFn`
in the SUNNONLINEARSOLVER API have been updated to remove unused input
parameters.

Added a new ``SUNNonlinearSolver`` implementation, ``SUNNonlinsol_PetscSNES``,
which interfaces to the PETSc SNES nonlinear solver API.

Added new Fortran 2003 interfaces for most SUNNONLINEARSOLVER modules. See the
:numref:`SUNDIALS.Fortran` section for more details.

**ARKODE changes**

The MRIStep module has been updated to support explicit, implicit, or ImEx
methods as the fast integrator using the ARKStep module. As a result some
function signatures have been changed including :c:func:`MRIStepCreate` which
now takes an ARKStep memory structure for the fast integration as an input.

Fixed a bug in the ARKStep time-stepping module that would result in an infinite
loop if the nonlinear solver failed to converge more than the maximum allowed times
during a single step.

Fixed a bug that would result in a "too much accuracy requested" error when
using fixed time step sizes with explicit methods in some cases.

Fixed a bug in ARKStep where the mass matrix linear solver setup function was
not called in the Matrix-free case.

Fixed a minor bug in ARKStep where an incorrect flag is reported when an
error occurs in the mass matrix setup or Jacobian-vector product setup
functions.

Fixed a memeory leak in FARKODE when not using the default nonlinear solver.

The reinitialization functions :c:func:`ERKStepReInit()`,
:c:func:`ARKStepReInit()`, and :c:func:`MRIStepReInit()` have been updated to
retain the minimum and maxiumum step size values from before reinitialization
rather than resetting them to the default values.

Removed extraneous calls to :c:func:`N_VMin()` for simulations where
the scalar valued absolute tolerance, or all entries of the
vector-valued absolute tolerance array, are strictly positive.  In
this scenario, ARKODE will remove at least one global reduction per
time step.

The ARKLS interface has been updated to only zero the Jacobian matrix before
calling a user-supplied Jacobian evaluation function when the attached linear
solver has type ``SUNLINEARSOLVER_DIRECT``.

A new linear solver interface function :c:func:`ARKLsLinSysFn` was added as an
alternative method for evaluating the linear system :math:`A = M - \gamma J`.

Added two new embedded ARK methods of orders 4 and 5 to ARKODE (from :cite:p:`KenCarp:19`).

Support for optional inequality constraints on individual components of the
solution vector has been added the ARKODE ERKStep and ARKStep modules. See
the descriptions of :c:func:`ERKStepSetConstraints()` and
:c:func:`ARKStepSetConstraints()` for more details. Note that enabling
constraint handling requires the NVECTOR operations :c:func:`N_VMinQuotient()`,
:c:func:`N_VConstrMask()`, and :c:func:`N_VCompare()` that were not previously
required by ARKODE.

Added two new 'Get' functions to ARKStep, :c:func:`ARKStepGetCurrentGamma()`,
and :c:func:`ARKStepGetCurrentState`, that may be useful to users who choose
to provide their own nonlinear solver implementation.

Add two new 'Set' functions to MRIStep, :c:func:`MRIStepSetPreInnerFn()` and
:c:func:`MRIStepSetPostInnerFn()` for performing communication or memory
transfers needed before or after the inner integration.

A new Fortran 2003 interface to ARKODE was added. This includes Fortran 2003 interfaces
to the ARKStep, ERKStep, and MRIStep time-stepping modules. See the
:numref:`SUNDIALS.Fortran` section for more details.



Changes in v3.1.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An additional NVECTOR implementation was added for the
Tpetra vector from the Trilinos library to facilitate interoperability
between SUNDIALS and Trilinos. This implementation is accompanied by
additions to user documentation and SUNDIALS examples.

A bug was fixed where a nonlinear solver object could be freed twice in some use
cases.

The ``EXAMPLES_ENABLE_RAJA`` CMake option has been removed. The option ``EXAMPLES_ENABLE_CUDA``
enables all examples that use CUDA including the RAJA examples with a CUDA back end
(if the RAJA NVECTOR is enabled).

The implementation header file `arkode_impl.h` is no longer installed. This means users
who are directly manipulating the ``ARKodeMem`` structure will need to update their code
to use ARKODE's public API.

Python is no longer required to run ``make test`` and ``make test_install``.

Fixed a bug in ``ARKodeButcherTable_Write`` when printing a Butcher table
without an embedding.

Changes in v3.0.2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Added information on how to contribute to SUNDIALS and a contributing agreement.

Changes in v3.0.1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A bug in ARKODE where single precision builds would fail to compile has been fixed.


Changes in v3.0.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ARKODE library has been entirely rewritten to support a modular
approach to one-step methods, which should allow rapid research and
development of novel integration methods without affecting existing
solver functionality.  To support this, the existing ARK-based methods
have been encapsulated inside the new ``ARKStep`` time-stepping
module. Two new time-stepping modules have been added:

* The ``ERKStep`` module provides an optimized implementation for explicit
  Runge--Kutta methods with reduced storage and number of calls to the ODE
  right-hand side function.

* The ``MRIStep`` module implements two-rate explicit-explicit multirate
  infinitesimal step methods utilizing different step sizes for slow
  and fast processes in an additive splitting.

This restructure has resulted in numerous small changes to the user
interface, particularly the suite of "Set" routines for user-provided
solver parameters and "Get" routines to access solver statistics,
that are now prefixed with the name of time-stepping module (e.g., ``ARKStep``
or ``ERKStep``) instead of ``ARKODE``.  Aside from affecting the names of these
routines, user-level changes have been kept to a minimum.  However, we recommend
that users consult both this documentation and the ARKODE example programs for
further details on the updated infrastructure.

As part of the ARKODE restructuring an :c:type:`ARKodeButcherTable` structure
has been added for storing Butcher tables. Functions for creating new Butcher
tables and checking their analytic order are provided along with other utility
routines. For more details see :numref:`ARKodeButcherTable`.

Two changes were made in the initial step size algorithm:

* Fixed an efficiency bug where an extra call to the right hand side function was made.

* Changed the behavior of the algorithm if the max-iterations case is hit.
  Before the algorithm would exit with the step size calculated on the
  penultimate iteration. Now it will exit with the step size calculated
  on the final iteration.

ARKODE's dense output infrastructure has been improved to support
higher-degree Hermite polynomial interpolants (up to degree 5) over
the last successful time step.

ARKODE's previous direct and iterative linear solver interfaces, ARKDLS and
ARKSPILS, have been merged into a single unified linear solver interface, ARKLS,
to support any valid SUNLINSOL module. This includes ``DIRECT`` and
``ITERATIVE`` types as well as the new ``MATRIX_ITERATIVE`` type. Details
regarding how ARKLS utilizes linear solvers of each type as well as discussion
regarding intended use cases for user-supplied SUNLinSol implementations are
included in the chapter :numref:`SUNLinSol`. All ARKODE examples programs and the
standalone linear solver examples have been updated to use the unified linear
solver interface.

The user interface for the new ARKLS module is very similar to the previous
ARKDLS and ARKSPILS interfaces. Additionally, we note that Fortran users will
need to enlarge their ``iout`` array of optional integer outputs, and update the
indices that they query for certain linear-solver-related statistics.

The names of all constructor routines for SUNDIALS-provided SUNLinSol
implementations have been updated to follow the naming convention
``SUNLinSol_*`` where ``*`` is the name of the linear solver. The new names are
``SUNLinSol_Band``, ``SUNLinSol_Dense``, ``SUNLinSol_KLU``,
``SUNLinSol_LapackBand``, ``SUNLinSol_LapackDense``, ``SUNLinSol_PCG``,
``SUNLinSol_SPBCGS``, ``SUNLinSol_SPFGMR``, ``SUNLinSol_SPGMR``,
``SUNLinSol_SPTFQMR``, and ``SUNLinSol_SuperLUMT``.  Solver-specific "set"
routine names have been similarly standardized.  To minimize challenges in user
migration to the new names, the previous routine names may still be used; these
will be deprecated in future releases, so we recommend that users migrate to the
new names soon. All ARKODE example programs and the standalone linear solver
examples have been updated to use the new naming convention.

The ``SUNBandMatrix`` constructor has been simplified to remove the
storage upper bandwidth argument.

SUNDIALS integrators have been updated to utilize generic nonlinear solver
modules defined through the SUNNONLINSOL API. This API will ease the addition of
new nonlinear solver options and allow for external or user-supplied nonlinear
solvers. The SUNNONLINSOL API and SUNDIALS provided modules are described in
:numref:`SUNNonlinSol` and follow the same object oriented design and
implementation used by the NVector, SUNMatrix, and SUNLinSol modules. Currently
two SUNNONLINSOL implementations are provided, SUNNonlinSol_Newton and
SUNNonlinSol_FixedPoint. These replicate the previous integrator specific
implementations of a Newton iteration and an accelerated fixed-point iteration,
respectively. Example programs using each of these nonlinear solver modules in a
standalone manner have been added and all ARKODE example programs have been
updated to use generic SUNNonlinSol modules.

As with previous versions, ARKODE will use the Newton solver (now
provided by SUNNonlinSol_Newton) by default.  Use of the
:c:func:`ARKStepSetLinear()` routine (previously named
``ARKodeSetLinear``) will indicate that the problem is
linearly-implicit, using only a single Newton iteration per implicit
stage.  Users wishing to switch to the accelerated fixed-point solver
are now required to create a SUNNonlinSol_FixedPoint object and attach
that to ARKODE, instead of calling the previous
``ARKodeSetFixedPoint`` routine.  See the documentation sections
:numref:`ARKODE.Usage.ARKStep.Skeleton`,
:numref:`ARKODE.Usage.ARKStep.NonlinearSolvers`, and
:numref:`SUNNonlinSol.FixedPoint` for further details, or the serial C
example program ``ark_brusselator_fp.c`` for an example.

Three fused vector operations and seven vector array operations have been added
to the NVECTOR API. These *optional* operations are disabled by default and may
be activated by calling vector specific routines after creating an NVector (see
:numref:`NVectors.Description` for more details). The new operations are intended
to increase data reuse in vector operations, reduce parallel communication on
distributed memory systems, and lower the number of kernel launches on systems
with accelerators. The fused operations are ``N_VLinearCombination``,
``N_VScaleAddMulti``, and ``N_VDotProdMulti``, and the vector array operations
are ``N_VLinearCombinationVectorArray``, ``N_VScaleVectorArray``,
``N_VConstVectorArray``, ``N_VWrmsNormVectorArray``,
``N_VWrmsNormMaskVectorArray``, ``N_VScaleAddMultiVectorArray``, and
``N_VLinearCombinationVectorArray``. If an NVector implementation defines any of
these operations as ``NULL``, then standard NVector operations will
automatically be called as necessary to complete the computation.

Multiple changes to the CUDA NVECTOR were made:

* Changed the ``N_VMake_Cuda`` function to take a host data pointer and a device
  data pointer instead of an ``N_VectorContent_Cuda`` object.

* Changed ``N_VGetLength_Cuda`` to return the global vector length instead of
  the local vector length.

* Added ``N_VGetLocalLength_Cuda`` to return the local vector length.

* Added ``N_VGetMPIComm_Cuda`` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``suncudavec``.

* Added the ability to set the ``cudaStream_t`` used for execution of the CUDA
  NVECTOR kernels. See the function ``N_VSetCudaStreams_Cuda``.

* Added ``N_VNewManaged_Cuda``, ``N_VMakeManaged_Cuda``, and ``N_VIsManagedMemory_Cuda``
  functions to accommodate using managed memory with the CUDA NVECTOR.

Multiple changes to the RAJA NVECTOR were made:

* Changed ``N_VGetLength_Raja`` to return the global vector length instead of
  the local vector length.

* Added ``N_VGetLocalLength_Raja`` to return the local vector length.

* Added ``N_VGetMPIComm_Raja`` to return the MPI communicator used.

* Removed the accessor functions in the namespace ``sunrajavec``.

A new NVECTOR implementation for leveraging OpenMP 4.5+ device offloading has
been added, NVECTOR_OpenMPDEV. See :numref:`NVectors.OpenMPDEV` for more details.


Changes in v2.2.1
^^^^^^^^^^^^^^^^^^^^^^^

Fixed a bug in the CUDA NVECTOR where the ``N_VInvTest`` operation could
write beyond the allocated vector data.

Fixed library installation path for multiarch systems. This fix changes the default
library installation path to ``CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR``
from ``CMAKE_INSTALL_PREFIX/lib``. ``CMAKE_INSTALL_LIBDIR`` is automatically
set, but is available as a CMAKE option that can modified.


Changes in v2.2.0
^^^^^^^^^^^^^^^^^^^^^^^

Fixed a problem with setting ``sunindextype`` which would occur with
some compilers (e.g. armclang) that did not define ``__STDC_VERSION__``.

Added hybrid MPI/CUDA and MPI/RAJA vectors to allow use of more than
one MPI rank when using a GPU system.  The vectors assume one GPU
device per MPI rank.

Changed the name of the RAJA NVECTOR library to
``libsundials_nveccudaraja.lib`` from
``libsundials_nvecraja.lib`` to better reflect that we only support CUDA
as a backend for RAJA currently.

Several changes were made to the build system:

* CMake 3.1.3 is now the minimum required CMake version.

* Deprecate the behavior of the ``SUNDIALS_INDEX_TYPE`` CMake option and
  added the ``SUNDIALS_INDEX_SIZE`` CMake option to select the ``sunindextype``
  integer size.

* The native CMake FindMPI module is now used to locate an MPI
  installation.

* If MPI is enabled and MPI compiler wrappers are not set, the build system
  will check if ``CMAKE_<language>_COMPILER`` can compile MPI programs before
  trying to locate and use an MPI installation.

* The previous options for setting MPI compiler wrappers and the executable
  for running MPI programs have been have been depreated. The new options that
  align with those used in native CMake FindMPI module are
  ``MPI_C_COMPILER``, ``MPI_CXX_COMPILER``, ``MPI_Fortran_COMPILER``,
  and ``MPIEXEC_EXECUTABLE``.

* When a Fortran name-mangling scheme is needed (e.g., ``ENABLE_LAPACK``
  is ``ON``) the build system will infer the scheme from the Fortran
  compiler. If a Fortran compiler is not available or the inferred or default
  scheme needs to be overridden, the advanced options
  ``SUNDIALS_F77_FUNC_CASE`` and ``SUNDIALS_F77_FUNC_UNDERSCORES`` can
  be used to manually set the name-mangling scheme and bypass trying to infer
  the scheme.

* Parts of the main CMakeLists.txt file were moved to new files in the
  ``src`` and ``example`` directories to make the CMake configuration file
  structure more modular.



Changes in v2.1.2
^^^^^^^^^^^^^^^^^^^^^^^

Updated the minimum required version of CMake to 2.8.12 and enabled
using rpath by default to locate shared libraries on OSX.

Fixed Windows specific problem where sunindextype was not correctly
defined when using 64-bit integers for the SUNDIALS index type. On Windows
sunindextype is now defined as the MSVC basic type ``__int64``.

Added sparse SUNMatrix "Reallocate" routine to allow specification of
the nonzero storage.

Updated the KLU SUNLinearSolver module to set constants for the two
reinitialization types, and fixed a bug in the full reinitialization
approach where the sparse SUNMatrix pointer would go out of scope on
some architectures.

Updated the "ScaleAdd" and "ScaleAddI" implementations in the
sparse SUNMatrix module to more optimally handle the case where the
target matrix contained sufficient storage for the sum, but had the
wrong sparsity pattern.  The sum now occurs in-place, by performing
the sum backwards in the existing storage.  However, it is still more
efficient if the user-supplied Jacobian routine allocates storage for
the sum :math:`I+\gamma J` or :math:`M+\gamma J` manually (with zero
entries if needed).

Changed LICENSE install path to ``instdir/include/sundials``.



Changes in v2.1.1
^^^^^^^^^^^^^^^^^^^

Fixed a potential memory leak in the SPGMR and SPFGMR linear solvers:
if "Initialize" was called multiple times then the solver memory was
reallocated (without being freed).

Fixed a minor bug in the ARKReInit routine, where a flag was
incorrectly set to indicate that the problem had been resized (instead
of just re-initialized).

Fixed C++11 compiler errors/warnings about incompatible use of string
literals.

Updated KLU SUNLinearSolver module to use a ``typedef`` for the
precision-specific solve function to be used (to avoid compiler
warnings).

Added missing typecasts for some ``(void*)`` pointers (again, to avoid
compiler warnings).

Bugfix in ``sunmatrix_sparse.c`` where we had used ``int`` instead of
``sunindextype`` in one location.

Added missing ``#include <stdio.h>`` in NVECTOR and SUNMATRIX header files.

Added missing prototype for ``ARKSpilsGetNumMTSetups``.

Fixed an indexing bug in the CUDA NVECTOR implementation of
``N_VWrmsNormMask`` and revised the RAJA NVECTOR implementation of
``N_VWrmsNormMask`` to work with mask arrays using values other than
zero or one. Replaced ``double`` with ``realtype`` in the RAJA vector
test functions.

Fixed compilation issue with GCC 7.3.0 and Fortran programs that do
not require a SUNMatrix or SUNLinearSolver module (e.g. iterative
linear solvers, explicit methods, fixed point solver, etc.).


Changes in v2.1.0
^^^^^^^^^^^^^^^^^^^

Added NVECTOR print functions that write vector data to a specified
file (e.g. ``N_VPrintFile_Serial``).

Added ``make test`` and ``make test_install`` options to the build
system for testing SUNDIALS after building with ``make`` and
installing with ``make install`` respectively.


Changes in v2.0.0
^^^^^^^^^^^^^^^^^^^

All interfaces to matrix structures and linear solvers have been
reworked, and all example programs have been updated.  The goal of the
redesign of these interfaces was to provide more encapsulation and
ease in interfacing custom linear solvers and interoperability with
linear solver libraries.

Specific changes include:

* Added generic SUNMATRIX module with three provided implementations:
  dense, banded and sparse.  These replicate previous SUNDIALS Dls and
  Sls matrix structures in a single object-oriented API.

* Added example problems demonstrating use of generic SUNMATRIX modules.

* Added generic SUNLINEARSOLVER module with eleven provided
  implementations: dense, banded, LAPACK dense, LAPACK band, KLU,
  SuperLU_MT, SPGMR, SPBCGS, SPTFQMR, SPFGMR, PCG.  These replicate
  previous SUNDIALS generic linear solvers in a single object-oriented
  API.

* Added example problems demonstrating use of generic SUNLINEARSOLVER modules.

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

Two additional NVECTOR implementations were added -- one for CUDA and
one for RAJA vectors.  These vectors are supplied to provide very
basic support for running on GPU architectures.  Users are advised
that these vectors both move all data to the GPU device upon
construction, and speedup will only be realized if the user also
conducts the right-hand-side function evaluation on the device. In
addition, these vectors assume the problem fits on one GPU. Further
information about RAJA, users are referred to the web site,
`https://software.llnl.gov/RAJA/ <https://software.llnl.gov/RAJA/>`_.
These additions are accompanied by additions to various interface
functions and to user documentation.

All indices for data structures were updated to a new ``sunindextype``
that can be configured to be a 32- or 64-bit integer data index type.
``sunindextype`` is defined to be ``int32_t`` or ``int64_t`` when
portable types are supported, otherwise it is defined as ``int`` or
``long int``. The Fortran interfaces continue to use ``long int`` for
indices, except for their sparse matrix interface that now uses the
new ``sunindextype``.  This new flexible capability for index types
includes interfaces to PETSc, *hypre*, SuperLU_MT, and KLU with either
32-bit or 64-bit capabilities depending how the user configures
SUNDIALS.

To avoid potential namespace conflicts, the macros defining
``booleantype`` values ``TRUE`` and ``FALSE`` have been changed to
``SUNTRUE`` and ``SUNFALSE`` respectively.

Temporary vectors were removed from preconditioner setup and solve
routines for all packages.  It is assumed that all necessary data
for user-provided preconditioner operations will be allocated and
stored in user-provided data structures.

The file ``include/sundials_fconfig.h`` was added.  This file contains
SUNDIALS type information for use in Fortran programs.

Added functions SUNDIALSGetVersion and SUNDIALSGetVersionNumber to get
SUNDIALS release version information at runtime.

The build system was expanded to support many of the xSDK-compliant keys.
The xSDK is a movement in scientific software to provide a foundation for the
rapid and efficient production of high-quality,
sustainable extreme-scale scientific applications.  More information can
be found at, `https://xsdk.info <https://xsdk.info>`_.

In addition, numerous changes were made to the build system.
These include the addition of separate ``BLAS_ENABLE`` and ``BLAS_LIBRARIES``
CMake variables, additional error checking during CMake configuration,
minor bug fixes, and renaming CMake options to enable/disable examples
for greater clarity and an added option to enable/disable Fortran 77 examples.
These changes included changing ``ENABLE_EXAMPLES`` to ``ENABLE_EXAMPLES_C``,
changing ``CXX_ENABLE`` to ``EXAMPLES_ENABLE_CXX``, changing ``F90_ENABLE`` to
``EXAMPLES_ENABLE_F90``, and adding an ``EXAMPLES_ENABLE_F77`` option.

Corrections and additions were made to the examples, to
installation-related files, and to the user documentation.




Changes in v1.1.0
^^^^^^^^^^^^^^^^^^^

We have included numerous bugfixes and enhancements since the
v1.0.2 release.

The bugfixes include:

* For each linear solver, the various solver performance counters are
  now initialized to 0 in both the solver specification function and
  in the solver's ``linit`` function.  This ensures that these solver
  counters are initialized upon linear solver instantiation as well as
  at the beginning of the problem solution.

* The choice of the method vs embedding the Billington and TRBDF2
  explicit Runge--Kutta methods were swapped, since in those the
  lower-order coefficients result in an A-stable method, while the
  higher-order coefficients do not.  This change results in
  significantly improved robustness when using those methods.

* A bug was fixed for the situation where a user supplies a vector of
  absolute tolerances, and also uses the vector Resize() functionality.

* A bug was fixed wherein a user-supplied Butcher table without an
  embedding is supplied, and the user is running with either fixed
  time steps (or they do adaptivity manually); previously this had
  resulted in an error since the embedding order was below 1.

* Numerous aspects of the documentation were fixed and/or clarified.


The feature changes/enhancements include:

* Two additional NVECTOR implementations were added -- one for Hypre
  (parallel) ParVector vectors, and one for PETSc vectors.  These
  additions are accompanied by additions to various interface
  functions and to user documentation.

* Each NVECTOR module now includes a function, ``N_VGetVectorID``,
  that returns the NVECTOR module name.

* A memory leak was fixed in the banded preconditioner and
  banded-block-diagonal preconditioner interfaces.  In addition,
  updates were done to return integers from linear solver and
  preconditioner 'free' routines.

* The Krylov linear solver Bi-CGstab was enhanced by removing a
  redundant dot product.  Various additions and corrections were made
  to the interfaces to the sparse solvers KLU and SuperLU_MT,
  including support for CSR format when using KLU.

* The ARKODE implicit predictor algorithms were updated: methods 2 and
  3 were improved slightly, a new predictor approach was added, and
  the default choice was modified.

* The underlying sparse matrix structure was enhanced to allow both
  CSR and CSC matrices, with CSR supported by the KLU linear solver
  interface.  ARKODE interfaces to the KLU solver from both C and
  Fortran were updated to enable selection of sparse matrix type, and a
  Fortran-90 CSR example program was added.

* The missing :c:func:`ARKSpilsGetNumMtimesEvals()` function was added
  -- this had been included in the previous documentation but had not
  been implemented.

* The handling of integer codes for specifying built-in ARKODE Butcher
  tables was enhanced.  While a global numbering system is still used,
  methods now have #defined names to simplify the user interface and to
  streamline incorporation of new Butcher tables into ARKODE.

* The maximum number of Butcher table stages was increased from 8 to
  15 to accommodate very high order methods, and an 8th-order adaptive
  ERK method was added.

* Support was added for the explicit and implicit methods in an
  additive Runge--Kutta method to utilize different stage times,
  solution and embedding coefficients, to support new SSP-ARK
  methods.

* The FARKODE interface was extended to include a routine to set
  scalar/array-valued residual tolerances, to support Fortran
  applications with non-identity mass-matrices.


Reading this User Guide
----------------------------

This user guide is a combination of general usage instructions and
specific example programs.  We expect that some readers will want to
concentrate on the general instructions, while others will refer
mostly to the examples, and the organization is intended to
accommodate both styles.

The structure of this document is as follows:

* In the next section we provide a thorough presentation of the
  underlying :ref:`mathematical algorithms <ARKODE.Mathematics>` used within
  the ARKODE family of solvers.

* We follow this with an overview of how the source code for both
  SUNDIALS and ARKODE are :ref:`organized <ARKODE.Organization>`.

* The largest section follows, providing a full account of how to use
  ARKODE's time-stepping modules, :ref:`ARKStep <ARKODE.Usage.ARKStep>`,
  :ref:`ERKStep <ARKODE.Usage.ERKStep>`, and :ref:`MRIStep <ARKODE.Usage.MRIStep>`,
  within C and C++ applications.  This section then includes additional
  information on how to use ARKODE from applications written in
  :ref:`Fortran <SUNDIALS.Fortran>`, as well as information on how to
  leverage :ref:`GPU accelerators within ARKODE <SUNDIALS.GPU>`.

* A much smaller section follows, describing ARKODE's
  :ref:`Butcher table structure <ARKodeButcherTable>`, that is used by
  both ARKStep and ERKStep.

* Subsequent sections discuss shared SUNDIALS features that are used
  by ARKODE:
  :ref:`vector data structures <NVectors>`,
  :ref:`matrix data structures <SUNMatrix>`,
  :ref:`linear solver data structures <SUNLinSol>`,
  :ref:`nonlinear solver data structures <SUNNonlinSol>`,
  :ref:`memory management utilities <SUNMemory>`,
  and the
  :ref:`installation procedure <Installation>`.

* The final sections catalog the full set of :ref:`ARKODE constants
  <ARKODE.Constants>`, that are used for both input specifications and return
  codes, and the full set of :ref:`Butcher tables <Butcher>` that are
  packaged with ARKODE.


SUNDIALS License and Notices
----------------------------

.. ifconfig:: package_name != 'super'

   .. include:: ../../../shared/LicenseReleaseNumbers.rst

.. ifconfig:: package_name == 'super'

   All SUNDIALS packages are released open source, under the BSD 3-Clause
   license for more details see the LICENSE and NOTICE files provided with all
   SUNDIALS packages.
