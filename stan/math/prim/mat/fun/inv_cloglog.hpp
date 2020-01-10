#ifndef STAN_MATH_PRIM_MAT_FUN_INV_CLOGLOG_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_CLOGLOG_HPP

#include <stan/math/prim/mat/fun/exp.hpp>
#include <stan/math/prim/fun/inv_cloglog.hpp>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv_cloglog() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 - exp(-exp(x)).
 */
struct inv_cloglog_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_cloglog(x);
  }
};

/**
 * Vectorized version of inv_cloglog().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto inv_cloglog(const T& x) {
  return apply_scalar_unary<inv_cloglog_fun, T>::apply(x);
}

/**
 * Version of inv_cloglog() that accepts Eigen Matrix or matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv_cloglog(const Eigen::MatrixBase<Derived>& x) {
  return (1 - exp(-exp(x.derived().array()))).matrix();
}

/**
 * Version of inv_cloglog() that accepts Eigen Array or array expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return 1 - exp(-exp()) applied to each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv_cloglog(const Eigen::ArrayBase<Derived>& x) {
  return 1 - exp(-exp(x.derived()));
}

}  // namespace math
}  // namespace stan

#endif
