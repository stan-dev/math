#ifndef STAN_MATH_PRIM_MAT_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_INV_SQRT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/inv_sqrt.hpp>

namespace stan {
namespace math {

/**
 * Structure to wrap inv_sqrt() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return 1 / sqrt of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_sqrt(x);
  }
};

/**
 * Vectorized version of inv_sqrt().
 *
 * @tparam T type of container
 * @param x container
 * @return 1 / sqrt of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto inv_sqrt(const T& x) {
  return apply_scalar_unary<inv_sqrt_fun, T>::apply(x);
}

/**
 * Version of inv_sqrt() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv_sqrt(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().rsqrt().matrix();
}

/**
 * Version of inv_sqrt() that accepts Eigen Array ar array expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv_sqrt(const Eigen::ArrayBase<Derived>& x) {
  return x.derived().rsqrt();
}

}  // namespace math
}  // namespace stan

#endif
