#ifndef STAN_MATH_PRIM_PROB_STD_NORMAL_LOG_QF_HPP
#define STAN_MATH_PRIM_PROB_STD_NORMAL_LOG_QF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inv_Phi.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_diff_exp.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse of the unit normal cumulative distribution function evaluated at
 * the log probability.
 *
 * @param log_p argument between -Inf and 0 inclusive
 * @return Real value of the inverse cdf for the standard normal distribution.
 */
inline double std_normal_log_qf(double log_p) {
  check_not_nan("std_normal_log_qf", "Log probability variable", log_p);
  check_less_or_equal("std_normal_log_qf", "Probability variable", log_p, 0);

  if (log_p == NEGATIVE_INFTY) {
    return NEGATIVE_INFTY;
  }
  if (log_p == 0) {
    return INFTY;
  }

  return internal::inv_Phi_impl<true>(log_p);
}

/**
 * Structure to wrap std_normal_log_qf() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable in range [-Inf, 0]
 * @return inverse of the unit normal CDF of x
 * @throw std::domain_error if x is not less than or equal to 0
 */
struct std_normal_log_qf_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return std_normal_log_qf(x);
  }
};

/**
 * A vectorized version of std_normal_log_qf() that accepts std::vectors, Eigen
 * Matrix/Array objects, or expressions, and containers of these.
 *
 * @tparam T type of container
 * @param x container
 * @return inverse unit normal CDF of each value in x
 * @throw std::domain_error if x is not less than or equal to 0
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto std_normal_log_qf(const T& x) {
  return apply_scalar_unary<std_normal_log_qf_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
