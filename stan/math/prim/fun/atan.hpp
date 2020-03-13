#ifndef STAN_MATH_PRIM_FUN_ATAN_HPP
#define STAN_MATH_PRIM_FUN_ATAN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap atan() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arctan of x in radians.
 */
struct atan_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::atan;
    return atan(x);
  }
};

/**
 * Vectorized version of atan().
 *
 * @tparam T type of container
 * @param x container
 * @return Arctan of each value in x, in radians.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline typename apply_scalar_unary<atan_fun, T>::return_t atan(const T& x) {
  return apply_scalar_unary<atan_fun, T>::apply(x);
}

/**
 * Version of atan() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Elementwise atan of members of container.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto atan(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().atan().matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
