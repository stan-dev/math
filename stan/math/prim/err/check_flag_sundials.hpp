#ifndef STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP

#include <kinsol/kinsol.h>
#include <cvodes/cvodes.h>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>

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
 * Throws an exception message when the function KINSol()
 * (call to the solver) fails. When the exception is caused
 * by a tuning parameter the user controls, gives a specific
 * error.
 *
 * @param flag Error flag
 * @throw <code>std::runtime_error</code> if the flag is negative.
 */
  inline void kinsol_check(int flag, const char* func_name) {
  std::ostringstream ss;
  if (flag < 0) {
    ss << func_name << " failed with error flag " << flag << ": "
       << KINGetReturnFlagName(flag) << ".";
    if (flag == -6) {
      throw std::domain_error(ss.str());
    } else {
      throw std::runtime_error(ss.str());
    }
  }
}

}  // namespace math
}  // namespace stan
#endif
