#ifndef STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FLAG_SUNDIALS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>

namespace stan {
namespace math {

/**
 * Throws a std::runtime_error exception when a Sundial function fails
 * (i.e. returns a negative flag)
 *
 * @param flag Error flag
 * @param func_name Name of the function that returned the flag
 * @throw <code>std::runtime_error</code> if the flag is negative
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
 *
 * @param flag Error flag
 * @param max_num_steps Maximum number of iterations the algebra solver
 *   should take before throwing an error
 * @throw <code>std::domain_error</code> if flag means maximum number of
 *   iterations exceeded in the algebra solver.
 * @throw <code>std::runtime_error</code> if the flag is negative for
 *   any other reason.
 */
inline void check_flag_kinsol(int flag,
                              long int max_num_steps) {  // NOLINT(runtime/int)
  std::ostringstream ss;
  if (flag == -6) {
    throw_domain_error("algebra_solver", "maximum number of iterations",
                       max_num_steps, "(", ") was exceeded in the solve.");
  } else if (flag < 0) {
    ss << "algebra_solver failed with error flag " << flag << ".";
    throw std::runtime_error(ss.str());
  }
}

}  // namespace math
}  // namespace stan
#endif
