#ifndef STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/domain_error.hpp>
#include <kinsol/kinsol.h>
#include <cvodes/cvodes.h>
#include <idas/idas.h>
#include <arkode/arkode.h>
#include <array>

namespace stan {
namespace math {

#define CHECK_CVODES_CALL(call) cvodes_check(call, #call)
#define CHECK_IDAS_CALL(call) idas_check(call, #call)
#define CHECK_KINSOL_CALL(call) kinsol_check(call, #call)
#define CHECK_ARKODE_CALL(call) arkode_check(call, #call)

/**
 * Map cvodes error flag to actually error msg. Messages follow the CVODES
 * documentation's table of return flags. Note that we
 * don't use CVGetReturnFlagName function to retrieve the constant
 * because sanitizer indicates it contains mem leak.
 *
 * @param flag
 *
 * @return error msg string constant and actually informative msg
 */
inline std::array<std::string, 2> cvodes_flag_msg(int flag) {
  switch (flag) {
    case CV_TOO_MUCH_WORK:
      return {"CV_TOO_MUCH_WORK",
              "The solver took mxstep internal steps but could not reach tout"};
    case CV_TOO_MUCH_ACC:
      return {"CV_TOO_MUCH_ACC",
              "The solver could not satisfy the accuracy demanded by the user "
              "for some internal step"};
    case CV_ERR_FAILURE:
      return {"CV_ERR_FAILURE",
              "Error test failures occurred too many times during one internal "
              "time step or minimum step size was reached"};
    case CV_CONV_FAILURE:
      return {"CV_CONV_FAILURE",
              "Convergence test failures occurred too many times during one "
              "internal time step or minimum step size was reached"};
    case CV_LINIT_FAIL:
      return {"CV_LINIT_FAIL",
              "The linear solver's initialization function failed"};
    case CV_LSETUP_FAIL:
      return {"CV_LSETUP_FAIL",
              "The linear solver's setup function failed in an "
              "unrecoverable manner"};
    case CV_LSOLVE_FAIL:
      return {"CV_LSOLVE_FAIL",
              "The linear solver's solve function failed in an "
              "unrecoverable manner"};
    case CV_RHSFUNC_FAIL:
      return {"CV_RHSFUNC_FAIL",
              "The right-hand side function failed in an unrecoverable manner"};
    case CV_FIRST_RHSFUNC_ERR:
      return {"CV_FIRST_RHSFUNC_ERR",
              "The right-hand side function failed at the first call"};
    case CV_REPTD_RHSFUNC_ERR:
      return {"CV_REPTD_RHSFUNC_ERR",
              "The right-hand side function had repeated recoverable errors"};
    case CV_UNREC_RHSFUNC_ERR:
      return {"CV_UNREC_RHSFUNC_ERR",
              "The right-hand side function had a recoverable error, but no "
              "recovery is possible"};
    case CV_RTFUNC_FAIL:
      return {"CV_RTFUNC_FAIL",
              "The rootfinding function failed in an unrecoverable manner"};
    case CV_NLS_INIT_FAIL:
      return {"CV_NLS_INIT_FAIL", "The nonlinear solver's init routine failed"};
    case CV_NLS_SETUP_FAIL:
      return {"CV_NLS_SETUP_FAIL",
              "The nonlinear solver's setup routine failed"};
    case CV_CONSTR_FAIL:
      return {"CV_CONSTR_FAIL",
              "The inequality constraints were violated and the solver was "
              "unable to recover"};
    case CV_MEM_FAIL:
      return {"CV_MEM_FAIL", "A memory allocation failed"};
    case CV_MEM_NULL:
      return {"CV_MEM_NULL", "The cvode_mem argument was NULL"};
    case CV_ILL_INPUT:
      return {"CV_ILL_INPUT", "One of the function inputs is illegal"};
    case CV_NO_MALLOC:
      return {"CV_NO_MALLOC",
              "The CVODE memory block was not allocated by a call to "
              "CVodeMalloc"};
    case CV_BAD_K:
      return {"CV_BAD_K",
              "The derivative order k is larger than the order used"};
    case CV_BAD_T:
      return {"CV_BAD_T", "The time t is outside the last step taken"};
    case CV_BAD_DKY:
      return {"CV_BAD_DKY", "The output derivative vector is NULL"};
    case CV_TOO_CLOSE:
      return {"CV_TOO_CLOSE",
              "The output and initial times are too close to each other"};
    case CV_NO_QUAD:
      return {"CV_NO_QUAD", "Quadrature integration was not activated"};
    case CV_QRHSFUNC_FAIL:
      return {"CV_QRHSFUNC_FAIL",
              "The quadrature right-hand side function failed in an "
              "unrecoverable manner"};
    case CV_FIRST_QRHSFUNC_ERR:
      return {"CV_FIRST_QRHSFUNC_ERR",
              "The quadrature right-hand side function failed at the first "
              "call"};
    case CV_REPTD_QRHSFUNC_ERR:
      return {"CV_REPTD_QRHSFUNC_ERR",
              "The quadrature right-hand side function had repeated "
              "recoverable errors"};
    case CV_UNREC_QRHSFUNC_ERR:
      return {"CV_UNREC_QRHSFUNC_ERR",
              "The quadrature right-hand side function had a recoverable "
              "error, but no recovery is possible"};
    case CV_NO_SENS:
      return {"CV_NO_SENS",
              "Forward sensitivity integration was not activated"};
    case CV_SRHSFUNC_FAIL:
      return {"CV_SRHSFUNC_FAIL",
              "The sensitivity right-hand side function failed in an "
              "unrecoverable manner"};
    case CV_FIRST_SRHSFUNC_ERR:
      return {"CV_FIRST_SRHSFUNC_ERR",
              "The sensitivity right-hand side function failed at the first "
              "call"};
    case CV_REPTD_SRHSFUNC_ERR:
      return {"CV_REPTD_SRHSFUNC_ERR",
              "The sensitivity right-hand side function had repeated "
              "recoverable errors"};
    case CV_UNREC_SRHSFUNC_ERR:
      return {"CV_UNREC_SRHSFUNC_ERR",
              "The sensitivity right-hand side function had a recoverable "
              "error, but no recovery is possible"};
    case CV_BAD_IS:
      return {"CV_BAD_IS",
              "The sensitivity index is larger than the number of "
              "sensitivities computed"};
    case CV_NO_QUADSENS:
      return {"CV_NO_QUADSENS",
              "Forward quadrature sensitivity integration was not activated"};
    case CV_QSRHSFUNC_FAIL:
      return {"CV_QSRHSFUNC_FAIL",
              "The quadrature sensitivity right-hand side function failed in "
              "an unrecoverable manner"};
    case CV_FIRST_QSRHSFUNC_ERR:
      return {"CV_FIRST_QSRHSFUNC_ERR",
              "The quadrature sensitivity right-hand side function failed at "
              "the first call"};
    case CV_REPTD_QSRHSFUNC_ERR:
      return {"CV_REPTD_QSRHSFUNC_ERR",
              "The quadrature sensitivity right-hand side function had "
              "repeated recoverable errors"};
    case CV_UNREC_QSRHSFUNC_ERR:
      return {"CV_UNREC_QSRHSFUNC_ERR",
              "The quadrature sensitivity right-hand side function had a "
              "recoverable error, but no recovery is possible"};
    case CV_CONTEXT_ERR:
      return {"CV_CONTEXT_ERR", "The SUNContext object is NULL"};
    case CV_NO_ADJ:
      return {"CV_NO_ADJ", "Adjoint module was not initialized"};
    case CV_NO_FWD:
      return {"CV_NO_FWD", "The forward integration was not yet performed"};
    case CV_NO_BCK:
      return {"CV_NO_BCK", "No backward problem was specified"};
    case CV_BAD_TB0:
      return {"CV_BAD_TB0",
              "The final time for the adjoint problem is outside the "
              "interval over which the forward problem was solved"};
    case CV_REIFWD_FAIL:
      return {"CV_REIFWD_FAIL",
              "Reinitialization of the forward problem failed at the first "
              "checkpoint"};
    case CV_FWD_FAIL:
      return {"CV_FWD_FAIL",
              "An error occurred during the integration of the forward "
              "problem"};
    case CV_GETY_BADT:
      return {"CV_GETY_BADT", "Wrong time in interpolation function"};
  }
  return {"CV_UNKNOWN_ERROR", "Unrecognized CVODES error flag"};
}

/**
 * Throws a std::domain_error exception when a Sundial function fails
 * (i.e. returns a negative flag)
 *
 * @param flag Error flag
 * @param func_name Name of the function that returned the flag
 * @throw <code>std::domain_error</code> if the flag is negative
 */
inline void cvodes_check(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag << ": \n"
       << cvodes_flag_msg(flag).at(1) << ".";
    throw std::domain_error(ss.str());
  }
}

/**
 * Map idas error flag to actually error msg. Messages follow the IDAS
 * documentation's table of return flags.
 *
 * @param flag
 *
 * @return error msg string constant and actually informative msg
 */
inline std::array<std::string, 2> idas_flag_msg(int flag) {
  switch (flag) {
    case IDA_TOO_MUCH_WORK:
      return {
          "IDA_TOO_MUCH_WORK",
          "The solver took mxstep internal steps but could not reach tout."};
    case IDA_TOO_MUCH_ACC:
      return {"IDA_TOO_MUCH_ACC",
              "The solver could not satisfy the accuracy demanded by the user "
              "for some internal step."};
    case IDA_ERR_FAIL:
      return {"IDA_ERR_FAIL",
              "Error test failures occurred too many times during one internal "
              "time step or minimum step size was reached."};
    case IDA_CONV_FAIL:
      return {"IDA_CONV_FAIL",
              "Convergence test failures occurred too many times during one "
              "internal time step or minimum step size was reached."};
    case IDA_LINIT_FAIL:
      return {"IDA_LINIT_FAIL",
              "The linear solver's initialization function failed."};
    case IDA_LSETUP_FAIL:
      return {"IDA_LSETUP_FAIL",
              "The linear solver's setup function failed in an unrecoverable "
              "manner."};
    case IDA_LSOLVE_FAIL:
      return {"IDA_LSOLVE_FAIL",
              "The linear solver's solve function failed in an unrecoverable "
              "manner."};
    case IDA_RES_FAIL:
      return {"IDA_RES_FAIL",
              "The user-provided residual function failed in an unrecoverable "
              "manner."};
    case IDA_REP_RES_ERR:
      return {"IDA_REP_RES_ERR",
              "The user-provided residual function repeatedly returned a "
              "recoverable error flag, but the solver was unable to recover."};
    case IDA_RTFUNC_FAIL:
      return {"IDA_RTFUNC_FAIL",
              "The rootfinding function failed in an unrecoverable manner."};
    case IDA_CONSTR_FAIL:
      return {"IDA_CONSTR_FAIL",
              "The inequality constraints were violated and the solver was "
              "unable to recover."};
    case IDA_FIRST_RES_FAIL:
      return {"IDA_FIRST_RES_FAIL",
              "The user-provided residual function failed recoverably on the "
              "first call."};
    case IDA_LINESEARCH_FAIL:
      return {"IDA_LINESEARCH_FAIL", "The line search failed."};
    case IDA_NO_RECOVERY:
      return {"IDA_NO_RECOVERY",
              "The residual function, linear solver setup function, or "
              "linear solver solve function had a recoverable failure, but "
              "IDACalcIC could not recover."};
    case IDA_NLS_INIT_FAIL:
      return {"IDA_NLS_INIT_FAIL",
              "The nonlinear solver's init routine failed."};
    case IDA_NLS_SETUP_FAIL:
      return {"IDA_NLS_SETUP_FAIL",
              "The nonlinear solver's setup routine failed."};
    case IDA_MEM_NULL:
      return {"IDA_MEM_NULL", "The ida_mem argument was NULL."};
    case IDA_MEM_FAIL:
      return {"IDA_MEM_FAIL", "A memory allocation failed."};
    case IDA_ILL_INPUT:
      return {"IDA_ILL_INPUT", "One of the function inputs is illegal."};
    case IDA_NO_MALLOC:
      return {"IDA_NO_MALLOC",
              "The IDAS memory was not allocated by a call to IDAInit."};
    case IDA_BAD_EWT:
      return {"IDA_BAD_EWT", "Zero value of some error weight component."};
    case IDA_BAD_K:
      return {"IDA_BAD_K", "The k-th derivative is not available."};
    case IDA_BAD_T:
      return {"IDA_BAD_T", "The time t is outside the last step taken."};
    case IDA_BAD_DKY:
      return {"IDA_BAD_DKY",
              "The vector argument where derivative should be stored is NULL."};
    case IDA_VECTOROP_ERR:
      return {"IDA_VECTOROP_ERR", "A vector operation failed."};
    case IDA_CONTEXT_ERR:
      return {"IDA_CONTEXT_ERR",
              "An error occurred in the SUNDIALS context object."};
    case IDA_NO_QUAD:
      return {"IDA_NO_QUAD", "Quadratures were not initialized."};
    case IDA_QRHS_FAIL:
      return {"IDA_QRHS_FAIL",
              "The user-provided right-hand side function for quadratures "
              "failed in an unrecoverable manner."};
    case IDA_FIRST_QRHS_ERR:
      return {"IDA_FIRST_QRHS_ERR",
              "The user-provided right-hand side function for quadratures "
              "failed in an unrecoverable manner on the first call."};
    case IDA_REP_QRHS_ERR:
      return {"IDA_REP_QRHS_ERR",
              "The user-provided right-hand side repeatedly returned a "
              "recoverable error flag, but the solver was unable to recover."};
    case IDA_NO_SENS:
      return {"IDA_NO_SENS", "Sensitivities were not initialized."};
    case IDA_SRES_FAIL:
      return {"IDA_SRES_FAIL",
              "The user-provided sensitivity residual function failed in an "
              "unrecoverable manner."};
    case IDA_REP_SRES_ERR:
      return {"IDA_REP_SRES_ERR",
              "The user-provided sensitivity residual function repeatedly "
              "returned a recoverable error flag, but the solver was "
              "unable to recover."};
    case IDA_BAD_IS:
      return {"IDA_BAD_IS", "The sensitivity identifier is not valid."};
    case IDA_NO_QUADSENS:
      return {"IDA_NO_QUADSENS",
              "Sensitivity-dependent quadratures were not initialized."};
    case IDA_QSRHS_FAIL:
      return {"IDA_QSRHS_FAIL",
              "The user-provided sensitivity-dependent quadrature right-hand "
              "side function failed in an unrecoverable manner."};
    case IDA_FIRST_QSRHS_ERR:
      return {"IDA_FIRST_QSRHS_ERR",
              "The user-provided sensitivity-dependent quadrature right-hand "
              "side function failed in an unrecoverable manner on the first "
              "call."};
    case IDA_REP_QSRHS_ERR:
      return {"IDA_REP_QSRHS_ERR",
              "The user-provided sensitivity-dependent quadrature right-hand "
              "side repeatedly returned a recoverable error flag, but the "
              "solver was unable to recover."};
    case IDA_NO_ADJ:
      return {"IDA_NO_ADJ",
              "The combined forward-backward problem has not been "
              "initialized."};
    case IDA_NO_FWD:
      return {"IDA_NO_FWD", "IDASolveF has not been previously called."};
    case IDA_NO_BCK:
      return {"IDA_NO_BCK", "No backward problem was specified."};
    case IDA_BAD_TB0:
      return {"IDA_BAD_TB0",
              "The desired output for backward problem is outside the "
              "interval over which the forward problem was solved."};
    case IDA_REIFWD_FAIL:
      return {"IDA_REIFWD_FAIL",
              "No checkpoint is available for this hot start."};
    case IDA_FWD_FAIL:
      return {"IDA_FWD_FAIL",
              "IDASolveB failed because IDASolve was unable to store data "
              "between two consecutive checkpoints."};
    case IDA_GETY_BADT:
      return {"IDA_GETY_BADT", "Wrong time in interpolation function."};
  }
  return {"IDA_UNKNOWN_ERROR", "Unrecognized IDAS error flag."};
}

inline void idas_check(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag << ": \n"
       << idas_flag_msg(flag).at(1);
    throw std::domain_error(ss.str());
  }
}

/**
 * Map ARKODE error flag to an error msg. Messages follow the ARKODE
 * documentation's table of shared output constants. The case labels use
 * the macros from arkode.h because the numeric values of several flags
 * differ between SUNDIALS releases.
 *
 * @param flag
 *
 * @return error msg string constant and actually informative msg
 */
inline std::array<std::string, 2> arkode_flag_msg(int flag) {
  switch (flag) {
    case ARK_TOO_MUCH_WORK:
      return {"ARK_TOO_MUCH_WORK",
              "The solver took mxstep internal steps but could not reach "
              "tout"};
    case ARK_TOO_MUCH_ACC:
      return {"ARK_TOO_MUCH_ACC",
              "The solver could not satisfy the accuracy demanded by the user "
              "for some internal step"};
    case ARK_ERR_FAILURE:
      return {"ARK_ERR_FAILURE",
              "Error test failures occurred too many times during one "
              "internal time step, or the minimum step size was reached"};
    case ARK_CONV_FAILURE:
      return {"ARK_CONV_FAILURE",
              "Convergence test failures occurred too many times during one "
              "internal time step, or the minimum step size was reached"};
    case ARK_LINIT_FAIL:
      return {"ARK_LINIT_FAIL",
              "The linear solver's initialization function failed"};
    case ARK_LSETUP_FAIL:
      return {"ARK_LSETUP_FAIL",
              "The linear solver's setup function failed in an "
              "unrecoverable manner"};
    case ARK_LSOLVE_FAIL:
      return {"ARK_LSOLVE_FAIL",
              "The linear solver's solve function failed in an "
              "unrecoverable manner"};
    case ARK_RHSFUNC_FAIL:
      return {"ARK_RHSFUNC_FAIL",
              "The right-hand side function failed in an unrecoverable "
              "manner"};
    case ARK_FIRST_RHSFUNC_ERR:
      return {"ARK_FIRST_RHSFUNC_ERR",
              "The right-hand side function failed at the first call"};
    case ARK_REPTD_RHSFUNC_ERR:
      return {"ARK_REPTD_RHSFUNC_ERR",
              "The right-hand side function had repeated recoverable errors"};
    case ARK_UNREC_RHSFUNC_ERR:
      return {"ARK_UNREC_RHSFUNC_ERR",
              "The right-hand side function had a recoverable error, but no "
              "recovery is possible"};
    case ARK_RTFUNC_FAIL:
      return {"ARK_RTFUNC_FAIL",
              "The rootfinding function failed in an unrecoverable manner"};
    case ARK_LFREE_FAIL:
      return {"ARK_LFREE_FAIL",
              "The linear solver's memory deallocation function failed"};
    case ARK_MASSINIT_FAIL:
      return {"ARK_MASSINIT_FAIL",
              "The mass matrix linear solver's initialization function "
              "failed"};
    case ARK_MASSSETUP_FAIL:
      return {"ARK_MASSSETUP_FAIL",
              "The mass matrix linear solver's setup function failed in an "
              "unrecoverable manner"};
    case ARK_MASSSOLVE_FAIL:
      return {"ARK_MASSSOLVE_FAIL",
              "The mass matrix linear solver's solve function failed in an "
              "unrecoverable manner"};
    case ARK_MASSFREE_FAIL:
      return {"ARK_MASSFREE_FAIL",
              "The mass matrix linear solver's memory deallocation "
              "function failed"};
    case ARK_MASSMULT_FAIL:
      return {"ARK_MASSMULT_FAIL",
              "The mass matrix-vector product function failed"};
    case ARK_CONSTR_FAIL:
      return {"ARK_CONSTR_FAIL",
              "The inequality constraint test failed repeatedly or failed "
              "with the minimum step size"};
    case ARK_MEM_FAIL:
      return {"ARK_MEM_FAIL", "A memory allocation failed"};
    case ARK_MEM_NULL:
      return {"ARK_MEM_NULL", "The arkode_mem argument was NULL"};
    case ARK_ILL_INPUT:
      return {"ARK_ILL_INPUT", "One of the function inputs is illegal"};
    case ARK_NO_MALLOC:
      return {"ARK_NO_MALLOC",
              "The ARKODE memory block was not allocated by a call to "
              "ARKStepCreate, ERKStepCreate, or MRIStepCreate"};
    case ARK_BAD_K:
      return {"ARK_BAD_K", "The derivative order k is larger than allowed"};
    case ARK_BAD_T:
      return {"ARK_BAD_T", "The time t is outside the last step taken"};
    case ARK_BAD_DKY:
      return {"ARK_BAD_DKY", "The output derivative vector is NULL"};
    case ARK_TOO_CLOSE:
      return {"ARK_TOO_CLOSE",
              "The output and initial times are too close to each other"};
    case ARK_VECTOROP_ERR:
      return {"ARK_VECTOROP_ERR",
              "An error occurred when calling an N_Vector routine"};
    case ARK_NLS_INIT_FAIL:
      return {"ARK_NLS_INIT_FAIL",
              "An error occurred when initializing a SUNNonlinSol module"};
    case ARK_NLS_SETUP_FAIL:
      return {"ARK_NLS_SETUP_FAIL",
              "A non-recoverable error occurred when setting up a "
              "SUNNonlinSol module"};
    case ARK_NLS_SETUP_RECVR:
      return {"ARK_NLS_SETUP_RECVR",
              "A recoverable error occurred when setting up a SUNNonlinSol "
              "module"};
    case ARK_NLS_OP_ERR:
      return {"ARK_NLS_OP_ERR",
              "An error occurred when calling a set/get routine in a "
              "SUNNonlinSol module"};
    case ARK_INNERSTEP_ATTACH_ERR:
      return {"ARK_INNERSTEP_ATTACH_ERR",
              "An error occurred when attaching the inner stepper module"};
    case ARK_INNERSTEP_FAIL:
      return {"ARK_INNERSTEP_FAIL",
              "An error occurred in the inner stepper module"};
    case ARK_OUTERTOINNER_FAIL:
      return {"ARK_OUTERTOINNER_FAIL",
              "An error occurred in the MRIStep pre inner integrator "
              "function"};
    case ARK_INNERTOOUTER_FAIL:
      return {"ARK_INNERTOOUTER_FAIL",
              "An error occurred in the MRIStep post inner integrator "
              "function"};
    case ARK_POSTPROCESS_STEP_FAIL:
      return {"ARK_POSTPROCESS_STEP_FAIL",
              "An error occurred in a step postprocessing function"};
    case ARK_POSTPROCESS_STAGE_FAIL:
      return {"ARK_POSTPROCESS_STAGE_FAIL",
              "An error occurred in a stage postprocessing function"};
    case ARK_USER_PREDICT_FAIL:
      return {"ARK_USER_PREDICT_FAIL",
              "An error occurred in a user predictor function"};
    case ARK_INTERP_FAIL:
      return {"ARK_INTERP_FAIL",
              "An error occurred in the ARKODE polynomial interpolation "
              "module"};
    case ARK_INVALID_TABLE:
      return {"ARK_INVALID_TABLE",
              "An invalid Butcher or MRI table was encountered"};
    case ARK_CONTEXT_ERR:
      return {"ARK_CONTEXT_ERR",
              "An error occurred with the SUNDIALS context object"};
    case ARK_UNRECOGNIZED_ERROR:
      return {"ARK_UNRECOGNIZED_ERROR", "An unknown error was encountered"};
  }
  return {"ARK_UNKNOWN_ERROR", "Unrecognized ARKODE error flag"};
}

inline void arkode_check(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag << ": \n"
       << arkode_flag_msg(flag).at(1) << ".";
    throw std::domain_error(ss.str());
  }
}

/**
 * Throws an exception message when the functions in KINSOL
 * fails. "KINGetReturnFlagName()" from SUNDIALS has a mem leak bug so
 * until it's fixed we cannot use it to extract flag error string.
 *
 * @param flag Error flag
 * @param func_name calling function name
 * @throw <code>std::runtime_error</code> if the flag is negative.
 */
inline void kinsol_check(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << "algebra_solver failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

/**
 * Throws an exception message when the KINSol() call fails.
 * When the exception is caused
 * by a tuning parameter the user controls, gives a specific
 * error.
 *
 * @param flag Error flag
 * @param func_name calling function name
 * @param max_num_steps max number of nonlinear iterations.
 * @throw <code>std::runtime_error</code> if the flag is negative.
 * @throw <code>std::domain_error</code> if the flag indicates max
 * number of steps is exceeded.
 */
inline void kinsol_check(int flag, const char* func_name,
                         long int max_num_steps) {  // NOLINT(runtime/int)
  std::ostringstream ss;
  if (flag == -6) {
    domain_error("algebra_solver", "maximum number of iterations",
                 max_num_steps, "(", ") was exceeded in the solve.");
  } else if (flag == -11) {
    ss << "The linear solver’s setup function failed "
       << "in an unrecoverable manner.";
    throw std::runtime_error(ss.str());
  } else if (flag < 0) {
    ss << "algebra_solver failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

}  // namespace math
}  // namespace stan
#endif
