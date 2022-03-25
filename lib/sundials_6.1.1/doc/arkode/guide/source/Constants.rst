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

.. _ARKODE.Constants:

===========================
Appendix: ARKODE Constants
===========================

Below we list all input and output constants used by the main solver,
timestepper, and linear solver modules, together with a short
description of their meaning.  :numref:`ARKODE.Constants.in_constants`
contains the ARKODE input constants, and :numref:`ARKODE.Constants.out_constants`
contains the ARKODE output constants.

.. _ARKODE.Constants.in_constants:
.. table:: ARKODE input constants
   :widths: 38 52

   +---------------------------------------------+-----------------------------------------------------------+
   | **Shared input constants**                                                                              |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_NORMAL`                         | Solver should return at a specified output time.          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_ONE_STEP`                       | Solver should return after each successful step.          |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **Full right-hand side evaluation constants**                                                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_FULLRHS_START`                  | Calling the full right-hand side function at the          |
   |                                             | start of the integration.                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_FULLRHS_END`                    | Calling the full right-hand side function at the end of   |
   |                                             | a step.                                                   |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_FULLRHS_OTHER`                  | Calling the full right-hand side function at the some     |
   |                                             | other point e.g., for dense output.                       |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **Interpolation module input constants**                                                                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_INTERP_HERMITE`                 | Specifies use of the Hermite polynomial interpolation     |
   |                                             | module (for non-stiff problems).                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_INTERP_LAGRANGE`                | Specifies use of the Lagrange polynomial interpolation    |
   |                                             | module (for stiff problems).                              |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARK_INTERP_MAX_DEGREE`              | Maximum possible interpolating polynomial degree.         |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **Explicit Butcher table specification**                                                                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_HEUN_EULER_2_1_2`            | Use the Heun-Euler-2-1-2 ERK method.                      |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_BOGACKI_SHAMPINE_4_2_3`      | Use the Bogacki-Shampine-4-2-3 ERK method.                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK324L2SA_ERK_4_2_3`        | Use the ARK-4-2-3 ERK method.                             |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ZONNEVELD_5_3_4`             | Use the Zonneveld-5-3-4 ERK method.                       |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK436L2SA_ERK_6_3_4`        | Use the ARK-6-3-4 ERK method.                             |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_SAYFY_ABURUB_6_3_4`          | Use the Sayfy-Aburub-6-3-4 ERK method.                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_CASH_KARP_6_4_5`             | Use the Cash-Karp-6-4-5 ERK method.                       |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_FEHLBERG_6_4_5`              | Use the Fehlberg-6-4-5 ERK method.                        |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_DORMAND_PRINCE_7_4_5`        | Use the Dormand-Prince-7-4-5 ERK method.                  |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK548L2SA_ERK_8_4_5`        | Use the ARK-8-4-5 ERK method.                             |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_VERNER_8_5_6`                | Use the Verner-8-5-6 ERK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_FEHLBERG_13_7_8`             | Use the Fehlberg-13-7-8 ERK method.                       |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_KNOTH_WOLKE_3_3`             | Use the Knoth-Wolke-3-3 ERK method.                       |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK437L2SA_ERK_7_3_4`        | Use the ARK-7-3-4 ERK method.                             |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK548L2SAb_ERK_8_4_5`       | Use the ARK-8-4-5b ERK method.                            |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_2`              | Use ARKStep's default second-order ERK method             |
   |                                             | (ARKODE_HEUN_EULER_2_1_2).                                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_3`              | Use ARKStep's default third-order ERK method              |
   |                                             | (ARKODE_BOGACKI_SHAMPINE_4_2_3).                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_4`              | Use ARKStep's default fourth-order ERK method             |
   |                                             | (ARKODE_ZONNEVELD_5_3_4).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_5`              | Use ARKStep's default fifth-order ERK method              |
   |                                             | (ARKODE_CASH_KARP_6_4_5).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_6`              | Use ARKStep's default sixth-order ERK method              |
   |                                             | (ARKODE_VERNER_8_5_6).                                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ERK_8`              | Use ARKStep's default eighth-order ERK method             |
   |                                             | (ARKODE_FEHLBERG_13_7_8).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_2`                  | Use ERKStep's default second-order ERK method             |
   |                                             | (ARKODE_HEUN_EULER_2_1_2).                                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_3`                  | Use ERKStep's default third-order ERK method              |
   |                                             | (ARKODE_BOGACKI_SHAMPINE_4_2_3).                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_4`                  | Use ERKStep's default fourth-order ERK method             |
   |                                             | (ARKODE_ZONNEVELD_5_3_4).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_5`                  | Use ERKStep's default fifth-order ERK method              |
   |                                             | (ARKODE_CASH_KARP_6_4_5).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_6`                  | Use ERKStep's default sixth-order ERK method              |
   |                                             | (ARKODE_VERNER_8_5_6).                                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ERKSTEP_DEFAULT_8`                  | Use ERKStep's default eighth-order ERK method             |
   |                                             | (ARKODE_FEHLBERG_13_7_8).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **Implicit Butcher table specification**                                                                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_SDIRK_2_1_2`                 | Use the SDIRK-2-1-2 SDIRK method.                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_BILLINGTON_3_3_2`            | Use the Billington-3-3-2 SDIRK method.                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_TRBDF2_3_3_2`                | Use the TRBDF2-3-3-2 ESDIRK method.                       |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_KVAERNO_4_2_3`               | Use the Kvaerno-4-2-3 ESDIRK method.                      |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK324L2SA_DIRK_4_2_3`       | Use the ARK-4-2-3 ESDIRK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_CASH_5_2_4`                  | Use the Cash-5-2-4 SDIRK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_CASH_5_3_4`                  | Use the Cash-5-3-4 SDIRK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_SDIRK_5_3_4`                 | Use the SDIRK-5-3-4 SDIRK method.                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_KVAERNO_5_3_4`               | Use the Kvaerno-5-3-4 ESDIRK method.                      |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK436L2SA_DIRK_6_3_4`       | Use the ARK-6-3-4 ESDIRK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_KVAERNO_7_4_5`               | Use the Kvaerno-7-4-5 ESDIRK method.                      |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK548L2SA_DIRK_8_4_5`       | Use the ARK-8-4-5 ESDIRK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK437L2SA_DIRK_7_3_4`       | Use the ARK-7-3-4 ESDIRK method.                          |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_ARK548L2SAb_DIRK_8_4_5`      | Use the ARK-8-4-5b ESDIRK method.                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_2`             | Use ARKStep's default second-order DIRK method            |
   |                                             | (ARKODE_SDIRK_2_1_2).                                     |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_3`             | Use ARKStep's default third-order DIRK method             |
   |                                             | (ARKODE_ARK324L2SA_DIRK_4_2_3).                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_4`             | Use ARKStep's default fourth-order DIRK method            |
   |                                             | (ARKODE_SDIRK_5_3_4).                                     |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_DIRK_5`             | Use ARKStep's default fifth-order DIRK method             |
   |                                             | (ARKODE_ARK548L2SA_DIRK_8_4_5).                           |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **ImEx Butcher table specification**                                                                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | ARKODE_ARK324L2SA_ERK_4_2_3 &               | Use the :index:`ARK-4-2-3 ARK method`.                    |
   | ARKODE_ARK324L2SA_DIRK_4_2_3                |                                                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | ARKODE_ARK436L2SA_ERK_6_3_4 &               | Use the :index:`ARK-6-3-4 ARK method`.                    |
   | ARKODE_ARK436L2SA_DIRK_6_3_4                |                                                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | ARKODE_ARK437L2SA_ERK_7_3_4 &               | Use the :index:`ARK-7-3-4 ARK method`.                    |
   | ARKODE_ARK437L2SA_DIRK_7_3_4                |                                                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | ARKODE_ARK548L2SA_ERK_8_4_5 &               | Use the :index:`ARK-8-4-5 ARK method`.                    |
   | ARKODE_ARK548L2SA_DIRK_8_4_5                |                                                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | ARKODE_ARK548L2SAb_ERK_8_4_5 &              | Use the :index:`ARK-8-4-5b ARK method`.                   |
   | ARKODE_ARK548L2SAb_DIRK_8_4_5               |                                                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_3` &     | Use ARKStep's default third-order ARK method              |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_3`       | (ARKODE_ARK324L2SA_ERK_4_2_3 and                          |
   |                                             | ARKODE_ARK324L2SA_DIRK_4_2_3).                            |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_4` &     | Use ARKStep's default fourth-order ARK method             |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_4`       | (ARKODE_ARK436L2SA_ERK_6_3_4 and                          |
   |                                             | ARKODE_ARK436L2SA_DIRK_6_3_4).                            |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKSTEP_DEFAULT_ARK_ETABLE_5` &     | Use ARKStep's default fifth-order ARK method              |
   | :index:`ARKSTEP_DEFAULT_ARK_ITABLE_5`       | (ARKODE_ARK548L2SA_ERK_8_4_5 and                          |
   |                                             | ARKODE_ARK548L2SA_DIRK_8_4_5).                            |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **MRI method types**                                                                                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_EXPLICIT`                   | Use an explicit (at the slow time scale) MRI method.      |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_IMPLICIT`                   | Use an implicit (at the slow time scale) MRI method.      |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_IMEX`                       | Use an ImEx (at the slow time scale) MRI method.          |
   +---------------------------------------------+-----------------------------------------------------------+
   |                                                                                                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | **MRI coupling table specification**                                                                    |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_MIS_MW3`                     | Use the Knoth-Wolke-3 MIS method.                         |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_MRI_GARK_ERK33a`             | Use the ERK33a MRI-GARK method.                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_MRI_GARK_ERK45a`             | Use the ERK45a MRI-GARK method.                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_MRI_GARK_IRK21a`             | Use the IRK21a MRI-GARK method.                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_MRI_GARK_ESDIRK34a`          | Use the ESDIRK34a MRI-GARK method.                        |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_MRI_GARK_ESDIRK46a`          | Use the ESDIRK46a MRI-GARK method.                        |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_IMEX_MRI_GARK3a`             | Use the IMEX-MRI-GARK3a method.                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_IMEX_MRI_GARK3b`             | Use the IMEX-MRI-GARK3b method.                           |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`ARKODE_IMEX_MRI_GARK4`              | Use the IMEX-MRI-GARK4 method.                            |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_TABLE_3`       | Use MRIStep's default 3rd-order explicit method           |
   |                                             | (MIS_MW3).                                                |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_EXPL_TABLE_4`       | Use MRIStep's default 4th-order explicit method           |
   |                                             | (MRI_GARK_ERK45a).                                        |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_TABLE_2`    | Use MRIStep's default 2nd-order solve-decoupled implicit  |
   |                                             | method (MRI_GARK_IRK21a).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_TABLE_3`    | Use MRIStep's default 3rd-order solve-decoupled implicit  |
   |                                             | method (MRI_GARK_ESDIRK34a).                              |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMPL_SD_TABLE_4`    | Use MRIStep's default 4th-order solve-decoupled implicit  |
   |                                             | method (MRI_GARK_ESDIRK46a).                              |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_TABLE_3`    | Use MRIStep's default 3rd-order solve-decoupled ImEx      |
   |                                             | method (IMEX_MRI_GARK3b).                                 |
   +---------------------------------------------+-----------------------------------------------------------+
   | :index:`MRISTEP_DEFAULT_IMEX_SD_TABLE_4`    | Use MRIStep's default 4th-order solve-decoupled ImEx      |
   |                                             | method (IMEX_MRI_GARK4).                                  |
   +---------------------------------------------+-----------------------------------------------------------+



.. _ARKODE.Constants.out_constants:
.. table:: ARKODE output constants
   :widths: 25 5 60

   +-------------------------------------+------+------------------------------------------------------------+
   | **Shared output constants**                                                                             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_SUCCESS`                | 0    | Successful function return.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TSTOP_RETURN`           | 1    | ARKODE succeeded by reaching the specified stopping point. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ROOT_RETURN`            | 2    | ARKODE succeeded and found one more more roots.            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_WARNING`                | 99   | ARKODE succeeded but an unusual situation occurred.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_MUCH_WORK`          | -1   | The solver took ``mxstep`` internal steps but could not    |
   |                                     |      | reach ``tout``.                                            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_MUCH_ACC`           | -2   | The solver could not satisfy the accuracy                  |
   |                                     |      | demanded by the user for some internal step.               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ERR_FAILURE`            | -3   | Error test failures occurred too many times during one     |
   |                                     |      | internal time step, or the minimum step size was reached.  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONV_FAILURE`           | -4   | Convergence test failures occurred too many times during   |
   |                                     |      | one internal time step, or the minimum step size was       |
   |                                     |      | reached.                                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LINIT_FAIL`             | -5   | The linear solver's initialization function failed.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LSETUP_FAIL`            | -6   | The linear solver's setup function failed in an            |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LSOLVE_FAIL`            | -7   | The linear solver's solve function failed in an            |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RHSFUNC_FAIL`           | -8   | The right-hand side function failed in an                  |
   |                                     |      | unrecoverable manner.                                      |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_FIRST_RHSFUNC_ERR`      | -9   | The right-hand side function failed at the first call.     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_REPTD_RHSFUNC_ERR`      | -10  | The right-hand side function had repeated recoverable      |
   |                                     |      | errors.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_UNREC_RHSFUNC_ERR`      | -11  | The right-hand side function had a recoverable error, but  |
   |                                     |      | no recovery is possible.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_RTFUNC_FAIL`            | -12  | The rootfinding function failed in an unrecoverable        |
   |                                     |      | manner.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_LFREE_FAIL`             | -13  | The linear solver's memory deallocation function failed.   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSINIT_FAIL`          | -14  | The mass matrix linear solver's initialization function    |
   |                                     |      | failed.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSSETUP_FAIL`         | -15  | The mass matrix linear solver's setup function failed in   |
   |                                     |      | an unrecoverable manner.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSSOLVE_FAIL`         | -16  | The mass matrix linear solver's solve function failed in   |
   |                                     |      | an unrecoverable manner.                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSFREE_FAIL`          | -17  | The mass matrix linear solver's memory deallocation        |
   |                                     |      | function failed.                                           |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MASSMULT_FAIL`          | -18  | The mass matrix-vector product function failed.            |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_CONSTR_FAIL`            | -19  | The inequality constraint test failed repeatedly or        |
   |                                     |      | failed with the minimum step size.                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MEM_FAIL`               | -20  | A memory allocation failed.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_MEM_NULL`               | -21  | The ``arkode_mem`` argument was ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_ILL_INPUT`              | -22  | One of the function inputs is illegal.                     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NO_MALLOC`              | -23  | The ARKODE memory block was not allocated by               |
   |                                     |      | a call to :c:func:`ARKStepCreate`,                         |
   |                                     |      | :c:func:`ERKStepCreate`, or :c:func:`MRIStepCreate`.       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_K`                  | -24  | The derivative order :math:`k` is larger than allowed.     |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_T`                  | -25  | The time :math:`t` is outside the last step taken.         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_BAD_DKY`                | -26  | The output derivative vector is ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_TOO_CLOSE`              | -27  | The output and initial times are too close to each other.  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_VECTOROP_ERR`           | -28  | An error occurred when calling an :c:type:`N_Vector`       |
   |                                     |      | routine.                                                   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_INIT_FAIL`          | -29  | An error occurred when initializing a SUNNonlinSol module. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_SETUP_FAIL`         | -30  | A non-recoverable error occurred when setting up a         |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_SETUP_RECVR`        | -31  | A recoverable error occurred when setting up a             |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_NLS_OP_ERR`             | -32  | An error occurred when calling a set/get routine in a      |
   |                                     |      | SUNNonlinSol module.                                       |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INNERSTEP_ATTACH_ERR`   | -33  | An error occurred when attaching the inner stepper module. |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INNERSTEP_FAIL`         | -34  | An error occurred in the inner stepper module.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_PREINNERFN_FAIL`        | -35  | An error occurred in the MRIStep pre inner integrator      |
   |                                     |      | function.                                                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_POSTINNERFN_FAIL`       | -36  | An error occurred in the MRIStep post inner integrator     |
   |                                     |      | function.                                                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INTERP_FAIL`            | -40  | An error occurred in the ARKODE polynomial interpolation   |
   |                                     |      | module.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_INVALID_TABLE`          | -41  | An invalid Butcher or MRI table was encountered.           |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARK_UNRECOGNIZED_ERROR`     | -99  | An unknown error was encountered.                          |
   +-------------------------------------+------+------------------------------------------------------------+
   |                                                                                                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | **ARKLS linear solver module output constants**                                                         |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUCCESS`              | 0    | Successful function return.                                |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MEM_NULL`             | -1   | The ``arkode_mem`` argument was ``NULL``.                  |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_LMEM_NULL`            | -2   | The ARKLS linear solver interface has not been             |
   |                                     |      | initialized.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_ILL_INPUT`            | -3   | The ARKLS solver interface is not compatible with          |
   |                                     |      | the current :c:type:`N_Vector` module, or an input value   |
   |                                     |      | was illegal.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MEM_FAIL`             | -4   | A memory allocation request failed.                        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_PMEM_NULL`            | -5   | The preconditioner module has not been initialized.        |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSMEM_NULL`         | -6   | The ARKLS mass-matrix linear solver interface has not been |
   |                                     |      | initialized.                                               |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_JACFUNC_UNRECVR`      | -7   | The Jacobian function failed in an unrecoverable manner.   |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_JACFUNC_RECVR`        | -8   | The Jacobian function had a recoverable error.             |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSFUNC_UNRECVR`     | -9   | The mass matrix function failed in an unrecoverable        |
   |                                     |      | manner.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_MASSFUNC_RECVR`       | -10  | The mass matrix function had a recoverable error.          |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUNMAT_FAIL`          | -11  | An error occurred with the current :c:type:`SUNMatrix`     |
   |                                     |      | module.                                                    |
   +-------------------------------------+------+------------------------------------------------------------+
   | :index:`ARKLS_SUNLS_FAIL`           | -12  | An error occurred with the current                         |
   |                                     |      | :c:type:`SUNLinearSolver` module.                          |
   +-------------------------------------+------+------------------------------------------------------------+



..
   Commented-out table rows:

      +-------------------------------------+------+------------------------------------------------------------+
      | :index:`ARK_POSTPROCESS_STEP_FAIL`  | -37  | An error occurred when calling the user-provided           |
      |                                     |      | step-based :c:func:`ARKPostProcessFn` routine.             |
      +-------------------------------------+------+------------------------------------------------------------+
      | :index:`ARK_POSTPROCESS_STAGE_FAIL` | -38  | An error occurred when calling the user-provided           |
      |                                     |      | stage-based :c:func:`ARKPostProcessFn` routine.            |
      +-------------------------------------+------+------------------------------------------------------------+
