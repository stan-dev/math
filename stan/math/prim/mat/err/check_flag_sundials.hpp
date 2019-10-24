#ifndef STAN_MATH_PRIM_MAT_ERR_CHECK_FLAG_SUNDIALS_HPP
#define STAN_MATH_PRIM_MAT_ERR_CHECK_FLAG_SUNDIALS_HPP

#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Throws an exception when a Sundial function fails
 * (i.e. returns a negative flag)
 */
inline void check_flag_sundials(int flag, const char* func_name) {
  if (flag < 0) {
    std::ostringstream ss;
    ss << func_name << " failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

/**
 * Throws an exception message when the function KINSol()
 * (call to the solver) fails. When the exception is caused
 * by a tuning parameter the user controls, gives a specific
 * error.
 */
inline void check_flag_kinsol(int flag,
                              long int max_num_steps) {  // NOLINT(runtime/int)
  std::ostringstream ss;
  if (flag == -6) {
    ss << "algebra_solver: max number of iterations: " << max_num_steps
       << " exceeded.";
    throw std::runtime_error(ss.str());
  } else if (flag < 0) {
    ss << "algebra_solver failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

}  // namespace math
}  // namespace stan
#endif
