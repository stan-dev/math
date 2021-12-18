.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2021, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _KINSOL.Constants:

****************
KINSOL Constants
****************

Below we list all input and output constants used by the main solver and linear
solver modules, together with their numerical values and a short description of
their meaning.

.. _KINSOL.Constants.kinsol_in_KINSOL.Constants:

KINSOL input constants
======================

.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|

.. table:: KINSOL Main Solver Input Constants

  +----------------------+--------+--------------------------------------+
  | Constant Name        | Value  | Description                          |
  +======================+========+======================================+
  | ``KIN_ETACHOICE1``   | 1      | Use Eisenstat and Walker Choice 1    |
  |                      |        | for :math:`\eta`.                    |
  +----------------------+--------+--------------------------------------+
  | ``KIN_ETACHOICE2``   | 2      | Use Eisenstat and Walker Choice 2    |
  |                      |        | for :math:`\eta`.                    |
  +----------------------+--------+--------------------------------------+
  | ``KIN_ETACONSTANT``  | 3      | Use constant value for :math:`\eta`. |
  +----------------------+--------+--------------------------------------+
  | ``KIN_NONE``         | 0      | Use Newton iteration.                |
  +----------------------+--------+--------------------------------------+
  | ``KIN_LINESEARCH``   | 1      | Use Newton iteration with linesearch |
  |                      |        | globalization.                       |
  +----------------------+--------+--------------------------------------+
  | ``KIN_PICARD``       | 2      | Use Picard iteration.                |
  +----------------------+--------+--------------------------------------+

.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|

.. table:: Iterative Linear Solver Constants

  +----------------------+--------+---------------------------------------+
  | Constant Name        | Value  | Description                           |
  +======================+========+=======================================+
  | ``SUN_PREC_NONE``    | 0      | No preconditioning                    |
  +----------------------+--------+---------------------------------------+
  | ``SUN_PREC_RIGHT``   | 2      | Preconditioning on the right.         |
  +----------------------+--------+---------------------------------------+
  | ``SUN_MODIFIED_GS``  | 1      | Use modified Gram-Schmidt procedure.  |
  +----------------------+--------+---------------------------------------+
  | ``SUN_CLASSICAL_GS`` | 2      | Use classical Gram-Schmidt procedure. |
  +----------------------+--------+---------------------------------------+

.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|

.. table:: Anderson Acceleration Orthogonalization Method Constants

  +---------------------+--------+---------------------------------------------+
  | Constant Name       | Value  | Description                                 |
  +=====================+========+=============================================+
  | ``KIN_ORTH_MGS``    | 0      | Use Modified Gram-Schmidt for Anderson      |
  |                     |        | acceleration.                               |
  +---------------------+--------+---------------------------------------------+
  | ``KIN_ORTH_ICWY``   | 1      | Use Inverse Compact WY Modified             |
  |                     |        | Gram-Schmidt for Anderson acceleration.     |
  +---------------------+--------+---------------------------------------------+
  | ``KIN_ORTH_CGS2``   | 2      | Use Classical Gram-Schmidt with             |
  |                     |        | Reorthogonalization (CGS-2) for Anderson    |
  |                     |        | Acceleration.                               |
  +---------------------+--------+---------------------------------------------+
  | ``KIN_ORTH_DCGS2``  | 3      | Use CGS-2 with Delayed Reorthogonalization  |
  |                     |        | for Anderson acceleration.                  |
  +---------------------+--------+---------------------------------------------+

.. _KINSOL.Constants.kinsol_out_KINSOL.Constants:

KINSOL output constants
=======================

.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|

.. table:: KINSOL Main Solver Output Constants

  +------------------------------+-------+-------------------------------------+
  | Constant Name                | Value | Description                         |
  +==============================+=======+=====================================+
  | ``KIN_SUCCESS``              | 0     | Successful function return.         |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_INITIAL_GUESS_OK``     | 1     | The initial user-supplied guess     |
  |                              |       | already satisfies the stopping      |
  |                              |       | criterion.                          |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_STEP_LT_STPTOL``       | 2     | The stopping tolerance on scaled    |
  |                              |       | step length was satisfied.          |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_WARNING``              | 99    | A non-fatal warning. The solver     |
  |                              |       | will continue.                      |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_MEM_NULL``             | -1    | The ``kin_mem`` argument was        |
  |                              |       | ``NULL``.                           |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_ILL_INPUT``            | -2    | One of the function inputs is       |
  |                              |       | illegal.                            |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_NO_MALLOC``            | -3    | The KINSOL memory was not allocated |
  |                              |       | by a call to ``KINMalloc``.         |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_MEM_FAIL``             | -4    | A memory allocation failed.         |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_LINESEARCH_NONCONV``   | -5    | The linesearch algorithm was unable |
  |                              |       | to find an iterate sufficiently     |
  |                              |       | distinct from the current iterate.  |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_MAXITER_REACHED``      | -6    | The maximum number of nonlinear     |
  |                              |       | iterations has been reached.        |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_MXNEWT_5X_EXCEEDED``   | -7    | Five consecutive steps have been    |
  |                              |       | taken that satisfy a scaled step    |
  |                              |       | length test.                        |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_LINESEARCH_BCFAIL``    | -8    | The linesearch algorithm was unable |
  |                              |       | to satisfy the                      |
  |                              |       | :math:`\beta`-condition for         |
  |                              |       | ``nbcfails`` iterations.            |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_LINSOLV_NO_RECOVERY``  | -9    | The user-supplied routine           |
  |                              |       | preconditioner slve function failed |
  |                              |       | recoverably, but the preconditioner |
  |                              |       | is already current.                 |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_LINIT_FAIL``           | -10   | The linear solver’s initialization  |
  |                              |       | function failed.                    |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_LSETUP_FAIL``          | -11   | The linear solver’s setup function  |
  |                              |       | failed in an unrecoverable manner.  |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_LSOLVE_FAIL``          | -12   | The linear solver’s solve function  |
  |                              |       | failed in an unrecoverable manner.  |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_SYSFUNC_FAIL``         | -13   | The system function failed in an    |
  |                              |       | unrecoverable manner.               |
  +------------------------------+-------+-------------------------------------+
  | ``KIN_FIRST_SYSFUNC_ERR``    | -14   | The system function failed with a   |
  |                              |       | recoverable error at the first call.|
  +------------------------------+-------+-------------------------------------+
  | ``KIN_REPTD_SYSFUNC_ERR``    | -15   | The system function had repeated    |
  |                              |       | recoverable errors.                 |
  +------------------------------+-------+-------------------------------------+


.. tabularcolumns:: |\Y{0.3}|\Y{0.1}|\Y{0.6}|

.. table:: KINLS Linear Solver Interface Output Constants

  +------------------------+--------+------------------------------------------+
  | Constant Name          | Value  | Description                              |
  +========================+========+==========================================+
  | ``KINLS_SUCCESS``      | 0      | Successful function return.              |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_MEM_NULL``     | -1     | The ``kin_mem`` argument was ``NULL``.   |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_LMEM_NULL``    | -2     | The KINLS linear solver has not been     |
  |                        |        | initialized.                             |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_ILL_INPUT``    | -3     | The KINLS solver is not compatible with  |
  |                        |        | the current ``N_Vector`` module, or an   |
  |                        |        | input value was illegal.                 |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_MEM_FAIL``     | -4     | A memory allocation request failed.      |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_PMEM_NULL``    | -5     | The preconditioner module has not been   |
  |                        |        | initialized.                             |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_JACFUNC_ERR``  | -6     | The Jacobian function failed             |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_SUNMAT_FAIL``  | -7     | An error occurred with the current       |
  |                        |        | ``SUNMatrix`` module.                    |
  +------------------------+--------+------------------------------------------+
  | ``KINLS_SUNLS_FAIL``   | -8     | An error occurred with the current       |
  |                        |        | ``SUNLinearSolver`` module.              |
  +------------------------+--------+------------------------------------------+
