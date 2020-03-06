#ifndef STAN_MATH_PRIM_FUN_INV_SQRT_HPP
#define STAN_MATH_PRIM_FUN_INV_SQRT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap 1 / sqrt(x) so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return inverse square root of x.
 */
struct inv_sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sqrt;
    return inv(sqrt(x));
  }
};

/**
 * Return the elementwise 1 / sqrt(x) of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam T type of container
 * @param x container
 * @return inverse square root of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto inv_sqrt(const T& x) {
  return apply_scalar_unary<inv_sqrt_fun, T>::apply(x);
}

/**
 * Version of inv_sqrt() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return inverse square root of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv_sqrt(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().rsqrt().matrix().eval();
}

/**
 * Version of inv_sqrt() that accepts Eigen Array or array expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return inverse square root of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto inv_sqrt(const Eigen::ArrayBase<Derived>& x) {
  return x.derived().rsqrt().eval();
}

}  // namespace math
}  // namespace stan

#endif
