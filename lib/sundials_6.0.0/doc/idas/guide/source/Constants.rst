.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _IDAS.Constants:

**************
IDAS Constants
**************

Below we list all input and output constants used by the main solver and
linear solver modules, together with their numerical values and a short
description of their meaning.

.. _IDAS.Constants.in_constants:

IDAS input constants
====================

.. table:: IDAS Input Constants
  :align: center

  +------------------------------------+-----+----------------------------------------------------+
  |    **IDAS main solver module**     |     |                                                    |
  +====================================+=====+====================================================+
  |                                    |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_NORMAL``                     | 1   | Solver returns at specified output time.           |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_ONE_STEP``                   | 2   | Solver returns after each successful step.         |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_SIMULTANEOUS``               | 1   | Simultaneous corrector forward sensitivity method. |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_STAGGERED``                  | 2   | Staggered corrector forward sensitivity method.    |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_CENTERED``                   | 1   | Central difference quotient approximation          |
  |                                    |     | (:math:`2^{nd}` order) of the sensitivity RHS.     |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_FORWARD``                    | 2   | Forward difference quotient approximation          |
  |                                    |     | (:math:`1^{st}` order) of the sensitivity RHS.     |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_YA_YDP_INIT``                | 1   | Compute :math:`y_a` and :math:`\dot{y}_d`, given   |
  |                                    |     | :math:`y_d`.                                       |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_Y_INIT``                     | 2   | Compute :math:`y`, given :math:`\dot{y}`.          |
  +------------------------------------+-----+----------------------------------------------------+
  |                                    |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  | **IDAS adjoint solver module**     |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  |                                    |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_HERMITE``                    | 1   | Use Hermite interpolation.                         |
  +------------------------------------+-----+----------------------------------------------------+
  | ``IDA_POLYNOMIAL``                 | 2   | Use variable-degree polynomial interpolation.      |
  +------------------------------------+-----+----------------------------------------------------+
  |                                    |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  | **Iterative linear solver module** |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  |                                    |     |                                                    |
  +------------------------------------+-----+----------------------------------------------------+
  | ``SUN_PREC_NONE``                  | 0   | No preconditioning                                 |
  +------------------------------------+-----+----------------------------------------------------+
  | ``SUN_PREC_LEFT``                  | 1   | Preconditioning on the left.                       |
  +------------------------------------+-----+----------------------------------------------------+
  | ``SUN_MODIFIED_GS``                | 1   | Use modified Gram-Schmidt procedure.               |
  +------------------------------------+-----+----------------------------------------------------+
  | ``SUN_CLASSICAL_GS``               | 2   | Use classical Gram-Schmidt procedure.              |
  +------------------------------------+-----+----------------------------------------------------+


.. _IDAS.Constants.out_constants:

IDAS output constants
=====================

.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|
.. table:: IDAS Output Constants
  :align: center

  +-----------------------------------+------+--------------------------------------------------------+
  |    **IDAS main solver module**    |      |                                                        |
  +===================================+======+========================================================+
  |                                   |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_SUCCESS``                   | 0    | Successful function return.                            |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_TSTOP_RETURN``              | 1    | ``IDASolve`` succeeded by reaching the specified       |
  |                                   |      | stopping point.                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_ROOT_RETURN``               | 2    | ``IDASolve`` succeeded and found one or more roots.    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_WARNING``                   | 99   | ``IDASolve`` succeeded but an unusual situation        |
  |                                   |      | occurred.                                              |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_TOO_MUCH_WORK``             | -1   | The solver took ``mxstep`` internal steps but could    |
  |                                   |      | not reach tout.                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_TOO_MUCH_ACC``              | -2   | The solver could not satisfy the accuracy demanded     |
  |                                   |      | by the user for some internal step.                    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_ERR_FAIL``                  | -3   | Error test failures occurred too many times during     |
  |                                   |      | one internal time step or minimum step size was        |
  |                                   |      | reached.                                               |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_CONV_FAIL``                 | -4   | Convergence test failures occurred too many times      |
  |                                   |      | during one internal time step or minimum step size     |
  |                                   |      | was reached.                                           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_LINIT_FAIL``                | -5   | The linear solver’s initialization function failed.    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_LSETUP_FAIL``               | -6   | The linear solver’s setup function failed in an        |
  |                                   |      | unrecoverable manner.                                  |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_LSOLVE_FAIL``               | -7   | The linear solver’s solve function failed in an        |
  |                                   |      | unrecoverable manner.                                  |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_RES_FAIL``                  | -8   | The user-provided residual function failed in an       |
  |                                   |      | unrecoverable manner.                                  |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_REP_RES_FAIL``              | -9   | The user-provided residual function repeatedly         |
  |                                   |      | returned a recoverable error flag, but the solver      |
  |                                   |      | was unable to recover.                                 |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_RTFUNC_FAIL``               | -10  | The rootfinding function failed in an unrecoverable    |
  |                                   |      | manner.                                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_CONSTR_FAIL``               | -11  | The inequality constraints were violated and the       |
  |                                   |      | solver was unable to recover.                          |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_FIRST_RES_FAIL``            | -12  | The user-provided residual function failed             |
  |                                   |      | recoverably on the first call.                         |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_LINESEARCH_FAIL``           | -13  | The line search failed.                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_RECOVERY``               | -14  | The residual function, linear solver setup function,   |
  |                                   |      | or linear solver solve function had a recoverable      |
  |                                   |      | failure, but ``IDACalcIC`` could not recover.          |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NLS_INIT_FAIL``             | -15  | The nonlinear solver’s init routine failed.            |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NLS_SETUP_FAIL``            | -16  | The nonlinear solver’s setup routine failed.           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_MEM_NULL``                  | -20  | The ``ida_mem`` argument was ``NULL``.                 |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_MEM_FAIL``                  | -21  | A memory allocation failed.                            |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_ILL_INPUT``                 | -22  | One of the function inputs is illegal.                 |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_MALLOC``                 | -23  | The IDAS memory was not allocated by a call to         |
  |                                   |      | ``IDAInit``.                                           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_BAD_EWT``                   | -24  | Zero value of some error weight component.             |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_BAD_K``                     | -25  | The :math:`k`-th derivative is not available.          |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_BAD_T``                     | -26  | The time :math:`t` is outside the last step taken.     |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_BAD_DKY``                   | -27  | The vector argument where derivative should be         |
  |                                   |      | stored is ``NULL``.                                    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_QUAD``                   | -30  | Quadratures were not initialized.                      |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_QRHS_FAIL``                 | -31  | The user-provided right-hand side function for         |
  |                                   |      | quadratures failed in an unrecoverable manner.         |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_FIRST_QRHS_ERR``            | -32  | The user-provided right-hand side function for         |
  |                                   |      | quadratures failed in an unrecoverable manner on the   |
  |                                   |      | first call.                                            |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_REP_QRHS_ERR``              | -33  | The user-provided right-hand side repeatedly           |
  |                                   |      | returned a recoverable error flag, but the solver      |
  |                                   |      | was unable to recover.                                 |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_SENS``                   | -40  | Sensitivities were not initialized.                    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_SRES_FAIL``                 | -41  | The user-provided sensitivity residual function        |
  |                                   |      | failed in an unrecoverable manner.                     |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_REP_SRES_ERR``              | -42  | The user-provided sensitivity residual function        |
  |                                   |      | repeatedly returned a recoverable error flag, but      |
  |                                   |      | the solver was unable to recover.                      |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_BAD_IS``                    | -43  | The sensitivity identifier is not valid.               |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_QUADSENS``               | -50  | Sensitivity-dependent quadratures were not             |
  |                                   |      | initialized.                                           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_QSRHS_FAIL``                | -51  | The user-provided sensitivity-dependent quadrature     |
  |                                   |      | right-hand side function failed in an unrecoverable    |
  |                                   |      | manner.                                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_FIRST_QSRHS_ERR``           | -52  | The user-provided sensitivity-dependent quadrature     |
  |                                   |      | right-hand side function failed in an unrecoverable    |
  |                                   |      | manner on the first call.                              |
  +-----------------------------------+------+--------------------------------------------------------+
  |                                   |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_REP_QSRHS_ERR``             | -53  | The user-provided sensitivity-dependent quadrature     |
  |                                   |      | right-hand side repeatedly returned a recoverable      |
  |                                   |      | error flag, but the solver was unable to recover.      |
  +-----------------------------------+------+--------------------------------------------------------+
  |                                   |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | **IDAS adjoint solver module**    |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  |                                   |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_ADJ``                    | -101 | The combined forward-backward problem has not been     |
  |                                   |      | initialized.                                           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_FWD``                    | -102 | ``IDASolveF`` has not been previously called.          |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_NO_BCK``                    | -103 | No backward problem was specified.                     |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_BAD_TB0``                   | -104 | The desired output for backward problem is outside     |
  |                                   |      | the interval over which the forward problem was        |
  |                                   |      | solved.                                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_REIFWD_FAIL``               | -105 | No checkpoint is available for this hot start.         |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_FWD_FAIL``                  | -106 | ``IDASolveB`` failed because ``IDASolve`` was unable   |
  |                                   |      | to store data between two consecutive checkpoints.     |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDA_GETY_BADT``                 | -107 | Wrong time in interpolation function.                  |
  +-----------------------------------+------+--------------------------------------------------------+
  |                                   |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | **IDALS linear solver interface** |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  |                                   |      |                                                        |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_SUCCESS``                 | 0    | Successful function return.                            |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_MEM_NULL``                | -1   | The ``ida_mem`` argument was ``NULL``.                 |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_LMEM_NULL``               | -2   | The IDALS linear solver has not been                   |
  |                                   |      | initialized.                                           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_ILL_INPUT``               | -3   | The IDALS solver is not compatible with the            |
  |                                   |      | current ``N_Vector`` module, or an input value was     |
  |                                   |      | illegal.                                               |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_MEM_FAIL``                | -4   | A memory allocation request failed.                    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_PMEM_NULL``               | -5   | The preconditioner module has not been initialized.    |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_JACFUNC_UNRECVR``         | -6   | The Jacobian function failed in an unrecoverable       |
  |                                   |      | manner.                                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_JACFUNC_RECVR``           | -7   | The Jacobian function had a recoverable error.         |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_SUNMAT_FAIL``             | -8   | An error occurred with the current ``SUNMatrix``       |
  |                                   |      | module.                                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_SUNLS_FAIL``              | -9   | An error occurred with the current ``SUNLinearSolver`` |
  |                                   |      | module.                                                |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_NO_ADJ``                  | -101 | The combined forward-backward problem has not been     |
  |                                   |      | initialized.                                           |
  +-----------------------------------+------+--------------------------------------------------------+
  | ``IDALS_LMEMB_NULL``              | -102 | The linear solver was not initialized for the          |
  |                                   |      | backward phase.                                        |
  +-----------------------------------+------+--------------------------------------------------------+
