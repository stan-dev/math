.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2002-2022, Lawrence Livermore National Security
   and Southern Methodist University.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _CVODE.Constants:

***************
CVODE Constants
***************

Below we list all input and output constants used by the main solver and
linear solver modules, together with their numerical values and a short
description of their meaning.  :numref:`CVODE.Constants.in_constants` contains
the CVODE input constants, and :numref:`CVODE.Constants.out_constants` contains
the CVODE output constants.


.. _CVODE.Constants.in_constants:
.. table:: CVODE input constants

   +----------------------+-----+-------------------------------------------------+
   | ``CV_ADAMS``         | 1   | Adams-Moulton linear multistep method.          |
   +----------------------+-----+-------------------------------------------------+
   | ``CV_BDF``           | 2   | BDF linear multistep method.                    |
   +----------------------+-----+-------------------------------------------------+
   | ``CV_NORMAL``        | 1   | Solver returns at specified output time.        |
   +----------------------+-----+-------------------------------------------------+
   | ``CV_ONE_STEP``      | 2   | Solver returns after each successful step.      |
   +----------------------+-----+-------------------------------------------------+
   | ``SUN_PREC_NONE``    | 0   | No preconditioning                              |
   +----------------------+-----+-------------------------------------------------+
   | ``SUN_PREC_LEFT``    | 1   | Preconditioning on the left only.               |
   +----------------------+-----+-------------------------------------------------+
   | ``SUN_PREC_RIGHT``   | 2   | Preconditioning on the right only.              |
   +----------------------+-----+-------------------------------------------------+
   | ``SUN_PREC_BOTH``    | 3   | Preconditioning on both the left and the right. |
   +----------------------+-----+-------------------------------------------------+
   | ``SUN_MODIFIED_GS``  | 1   | Use modified Gram-Schmidt procedure.            |
   +----------------------+-----+-------------------------------------------------+
   | ``SUN_CLASSICAL_GS`` | 2   | Use classical Gram-Schmidt procedure.           |
   +----------------------+-----+-------------------------------------------------+


.. _CVODE.Constants.out_constants:
.. table:: CVODE output constants
   :widths: 25 5 60

   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | **General outputs**                                                                                                       |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_SUCCESS``             | 0   | Successful function return.                                                            |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_TSTOP_RETURN``        | 1   | CVode succeeded by reaching the specified stopping point.                              |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_ROOT_RETURN``         | 2   | CVode succeeded and found one or more roots.                                           |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_WARNING``             | 99  | CVode succeeded but an unusual situation occurred.                                     |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_TOO_MUCH_WORK``       | -1  | The solver took ``mxstep`` internal steps but could not reach tout.                    |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_TOO_MUCH_ACC``        | -2  | The solver could not satisfy the accuracy demanded by the user for some internal step. |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_ERR_FAILURE``         | -3  | Error test failures occurred too many times during one                                 |
   |                            |     | internal time step or minimum step size was reached.                                   |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_CONV_FAILURE``        | -4  | Convergence test failures occurred too many times during one                           |
   |                            |     | internal time step or minimum step size was reached.                                   |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_LINIT_FAIL``          | -5  | The linear solver's initialization function failed.                                    |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_LSETUP_FAIL``         | -6  | The linear solver's setup function failed in an unrecoverable manner.                  |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_LSOLVE_FAIL``         | -7  | The linear solver's solve function failed in an unrecoverable manner.                  |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_RHSFUNC_FAIL``        | -8  | The right-hand side function failed in an unrecoverable manner.                        |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_FIRST_RHSFUNC_ERR``   | -9  | The right-hand side function failed at the first call.                                 |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_REPTD_RHSFUNC_ERR``   | -10 | The right-hand side function had repetead recoverable errors.                          |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_UNREC_RHSFUNC_ERR``   | -11 | The right-hand side function had a recoverable error, but no recovery is possible.     |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_RTFUNC_FAIL``         | -12 | The rootfinding function failed in an unrecoverable manner.                            |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_NLS_INIT_FAIL``       | -13 | The nonlinear solver's init routine failed.                                            |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_NLS_SETUP_FAIL``      | -14 | The nonlinear solver's setup routine failed.                                           |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_CONSTR_FAIL``         | -15 | The inequality constraints were violated and the solver was unable to recover.         |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_NLS_FAIL``            | -16 | The nonlinear solver failed in an unrecoverable manner.                                |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_MEM_FAIL``            | -20 | A memory allocation failed.                                                            |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_MEM_NULL``            | -21 | The ``cvode_mem`` argument was ``NULL``.                                               |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_ILL_INPUT``           | -22 | One of the function inputs is illegal.                                                 |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_NO_MALLOC``           | -23 | The CVODE memory block was not allocated by a call to ``CVodeMalloc``.                 |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_BAD_K``               | -24 | The derivative order $k$ is larger than the order used.                                |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_BAD_T``               | -25 | The time $t$ is outside the last step taken.                                           |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_BAD_DKY``             | -26 | The output derivative vector is ``NULL``.                                              |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_TOO_CLOSE``           | -27 | The output and initial times are too close to each other.                              |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_VECTOROP_ERR``        | -28 | A vector operation failed.                                                             |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_PROJ_MEM_NULL``       | -29 | The projection memory was ``NULL``.                                                    |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_PROJFUNC_FAIL``       | -30 | The projection function failed in an unrecoverable manner.                             |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CV_REPTD_PROJFUNC_ERR``  | -31 | The projection function had repeated recoverable errors.                               |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | **CVLS linear solver interface outputs**                                                                                  |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_SUCCESS``           | 0   | Successful function return.                                                            |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_MEM_NULL``          | -1  | The ``cvode_mem`` argument was ``NULL``.                                               |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_LMEM_NULL``         | -2  | The CVLS linear solver has not been initialized.                                       |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_ILL_INPUT``         | -3  | The CVLS solver is not compatible with the current ``N_Vector`` module.                |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_MEM_FAIL``          | -4  | A memory allocation request failed.                                                    |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_PMEM_NULL``         | -5  | The preconditioner module has not been initialized.                                    |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_JACFUNC_UNRECVR``   | -6  | The Jacobian function failed in an unrecoverable manner.                               |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_JACFUNC_RECVR``     | -7  | The Jacobian function had a recoverable error.                                         |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_SUNMAT_FAIL``       | -8  | An error occurred with the current {\sunmatrix} module.                                |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVLS_SUNLS_FAIL``        | -9  | An error occurred with the current {\sunlinsol} module.                                |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | **CVDIAG linear solver outputs**                                                                                          |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_SUCCESS``         | 0   | Successful function return.                                                            |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_MEM_NULL``        | -1  | The ``cvode_mem`` argument was ``NULL``.                                               |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_LMEM_NULL``       | -2  | The CVDIAG linear solver has not been initialized.                                     |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_ILL_INPUT``       | -3  | The CVDIAG solver is not compatible with the current ``N_Vector`` module.              |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_MEM_FAIL``        | -4  | A memory allocation request failed.                                                    |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_INV_FAIL``        | -5  | A diagonal element of the Jacobian was 0.                                              |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_RHSFUNC_UNRECVR`` | -6  | The right-hand side function failed in an unrecoverable manner.                        |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
   | ``CVDIAG_RHSFUNC_RECVR``   | -7  | The right-hand side function had a recoverable error.                                  |
   +----------------------------+-----+----------------------------------------------------------------------------------------+
