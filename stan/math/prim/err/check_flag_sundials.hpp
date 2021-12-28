#ifndef STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/domain_error.hpp>
#include <kinsol/kinsol.h>
#include <cvodes/cvodes.h>

namespace stan {
namespace math {

#define CHECK_CVODES_CALL(call) cvodes_check(call, #call)
#define CHECK_KINSOL_CALL(call) kinsol_check(call, #call)

/**
 * Map cvodes error flag to acutally error msg. The most frequent
 * errors are put at the top. An alternative would be to use std::map
 * but in our case the difference would be negligible. Note that we
 * don't use CVGetReturnFlagName function to retrieve the constant
 * because sanitizer indicates it contains mem leak.
 *
 * @param flag
 *
 * @return error msg string constant and actuall informative msg
 */
inline std::array<std::string, 2> cvodes_flag_msg(int flag) {
  std::array<std::string, 2> msg;
  switch (flag) {
    case -1:
      msg = {"CV_TOO_MUCH_WORK",
             "The solver took mxstep internal steps but could not reach tout"};
      break;  // NOLINT
    case -2:
      msg = {"CV_TOO_MUCH_ACC",
             "The solver could not satisfy the accuracy demanded by the user "
             "for some internal step"};
      break;  // NOLINT
    case -3:
      msg = {"CV_ERR_FAILURE",
             "Error test failures occurred too many times during one internal "
             "time step or minimum step size was reached"};
      break;  // NOLINT
    case -4:
      msg = {"CV_CONV_FAILURE",
             "Convergence test failures occurred too many times during one "
             "internal time step or minimum step size was reached"};
      break;  // NOLINT
    case -8:
      msg = {"CV_RHSFUNC_FAIL",
             "The right-hand side function failed in an unrecoverable manner"};
      break;  // NOLINT
    case -9:
      msg = {"CV_FIRST_RHSFUNC_ERR",
             "The right-hand side function failed at the first call"};
      break;  // NOLINT
    case -10:
      msg = {"CV_REPTD_RHSFUNC_ERR",
             "The right-hand side function had repetead recoverable errors"};
      break;  // NOLINT
    case -11:
      msg = {"CV_UNREC_RHSFUNC_ERR",
             "The right-hand side function had a recoverable error, but no "
             "recovery is possible"};
      break;  // NOLINT
    case -27:
      msg = {"CV_TOO_CLOSE",
             "The output and initial times are too close to each other"};
      break;  // NOLINT
    default:
      switch (flag) {
        case -5:
          msg = {"CV_LINIT_FAIL",
                 "The linear solver's initialization function failed"};
          break;  // NOLINT
        case -6:
          msg = {"CV_LSETUP_FAIL",
                 "The linear solver's setup function failed in an "
                 "unrecoverable manner"};
          break;  // NOLINT
        case -7:
          msg = {"CV_LSOLVE_FAIL",
                 "The linear solver's solve function failed in an "
                 "unrecoverable manner"};
          break;  // NOLINT
        case -20:
          msg = {"CV_MEM_FAIL", "A memory allocation failed"};
          break;  // NOLINT
        case -21:
          msg = {"CV_MEM_NULL", "The cvode_mem argument was NULL"};
          break;  // NOLINT
        case -22:
          msg = {"CV_ILL_INPUT", "One of the function inputs is illegal"};
          break;  // NOLINT
        case -23:
          msg = {"CV_NO_MALLOC",
                 "The CVODE memory block was not allocated by a call to "
                 "CVodeMalloc"};
          break;  // NOLINT
        case -24:
          msg = {"CV_BAD_K",
                 "The derivative order k is larger than the order used"};
          break;  // NOLINT
        case -25:
          msg = {"CV_BAD_T", "The time t s outside the last step taken"};
          break;  // NOLINT
        case -26:
          msg = {"CV_BAD_DKY", "The output derivative vector is NULL"};
          break;  // NOLINT
        case -40:
          msg = {"CV_BAD_IS",
                 "The sensitivity index is larger than the number of "
                 "sensitivities computed"};
          break;  // NOLINT
        case -41:
          msg = {"CV_NO_SENS",
                 "Forward sensitivity integration was not activated"};
          break;  // NOLINT
        case -42:
          msg = {"CV_SRHSFUNC_FAIL",
                 "The sensitivity right-hand side function failed in an "
                 "unrecoverable manner"};
          break;  // NOLINT
        case -43:
          msg = {"CV_FIRST_SRHSFUNC_ER",
                 "The sensitivity right-hand side function failed at the first "
                 "call"};
          break;  // NOLINT
        case -44:
          msg = {"CV_REPTD_SRHSFUNC_ER",
                 "The sensitivity ight-hand side function had repetead "
                 "recoverable errors"};
          break;  // NOLINT
        case -45:
          msg = {"CV_UNREC_SRHSFUNC_ER",
                 "The sensitivity right-hand side function had a recoverable "
                 "error, but no recovery is possible"};
          break;  // NOLINT
        case -101:
          msg = {"CV_ADJMEM_NULL", "The cvadj_mem argument was NULL"};
          break;  // NOLINT
        case -103:
          msg = {"CV_BAD_TB0",
                 "The final time for the adjoint problem is outside the "
                 "interval over which the forward problem was solved"};
          break;  // NOLINT
        case -104:
          msg = {"CV_BCKMEM_NULL",
                 "The cvodes memory for the backward problem was not created"};
          break;  // NOLINT
        case -105:
          msg = {"CV_REIFWD_FAIL",
                 "Reinitialization of the forward problem failed at the first "
                 "checkpoint"};
          break;  // NOLINT
        case -106:
          msg = {
              "CV_FWD_FAIL",
              "An error occured during the integration of the forward problem"};
          break;  // NOLINT
        case -107:
          msg = {"CV_BAD_ITASK", "Wrong task for backward integration"};
          break;  // NOLINT
        case -108:
          msg = {"CV_BAD_TBOUT",
                 "The desired output time is outside the interval over which "
                 "the forward problem was solved"};
          break;  // NOLINT
        case -109:
          msg = {"CV_GETY_BADT", "Wrong time in interpolation function"};
          break;  // NOLINT
      }
  }
  return msg;
}

/**
 * Throws a std::runtime_error exception when a Sundial function fails
 * (i.e. returns a negative flag)
 *
 * @param flag Error flag
 * @param func_name Name of the function that returned the flag
 * @throw <code>std::runtime_error</code> if the flag is negative
 */
inline void cvodes_check(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag << ": \n"
       << cvodes_flag_msg(flag).at(1) << ".";
    if (flag == -1 || flag == -4) {
      throw std::domain_error(ss.str());
    } else {
      throw std::runtime_error(ss.str());
    }
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
