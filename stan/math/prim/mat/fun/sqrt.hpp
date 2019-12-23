#ifndef STAN_MATH_PRIM_MAT_FUN_SQRT_HPP
#define STAN_MATH_PRIM_MAT_FUN_SQRT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sqrt() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Square root of x.
 */
struct sqrt_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sqrt;
    return sqrt(x);
  }
};

/**
 * Vectorized version of sqrt().
 *
 * @tparam T type of container
 * @param x container
 * @return Square root of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto sqrt(const T& x) {
  return apply_scalar_unary<sqrt_fun, T>::apply(x);
}

/**
 * Version of sqrt() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Square root of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto sqrt(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().sqrt().matrix();
}

}  // namespace math
}  // namespace stan

#endif
