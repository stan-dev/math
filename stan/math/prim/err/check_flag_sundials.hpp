#ifndef STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP

#include <kinsol/kinsol.h>
#include <cvodes/cvodes.h>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/domain_error.hpp>

namespace stan {
namespace math {

#define CHECK_CVODES_CALL(call) cvodes_check(call, #call)
#define CHECK_KINSOL_CALL(call) kinsol_check(call, #call)

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
    ss << func_name << " failed with error flag " << flag << ": "
       << CVodeGetReturnFlagName(flag) << ".";
    if (flag == -1 || flag == -4) {
      throw std::domain_error(ss.str());
    } else {
      throw std::runtime_error(ss.str());
    }
  }
}

/**
 * Throws an exception message when the functions in KINSOL
 * fails. When the exception is caused
 * by a tuning parameter the user controls, gives a specific
 * error. "KINGetReturnFlagName()" from sundials has a mem leak bug so
 * until it's fixed we cannot use it to extract flag error string.
 *
 * @param flag Error flag
 * @param func_name calling function name
 * @throw <code>std::runtime_error</code> if the flag is negative.
 */
inline void kinsol_check(int flag, const char* func_name) {
  std::ostringstream ss;
  if (flag < 0) {
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
 * @param max_num_steps max number of nonlinear iters
 * @throw <code>std::runtime_error</code> if the flag is negative.
 * @throw <code>std::domain_error</code> if the flag indicates max
 * number of steps is exceeded..
 */
inline void kinsol_check(int flag, const char* func_name,
                         long int max_num_steps) {  // NOLINT(runtime/int)
  std::ostringstream ss;
  if (flag == -6) {
    domain_error("algebra_solver", "maximum number of iterations",
                 max_num_steps, "(", ") was exceeded in the solve.");
  } else if (flag == -11) {
    ss << "The linear solver’s setup function failed in an unrecoverable "
          "manner.";
    throw std::runtime_error(ss.str());
  } else if (flag < 0) {
    ss << "algebra_solver failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

}  // namespace math
}  // namespace stan
#endif
