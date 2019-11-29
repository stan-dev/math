#ifndef STAN_MATH_PRIM_MAT_FUN_ASIN_HPP
#define STAN_MATH_PRIM_MAT_FUN_ASIN_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap asin() so it can be vectorized.
 * @param x Argument variable.
 * @tparam T Argument type.
 * @return Arcsine of x in radians.
 */
struct asin_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::asin;
    return asin(x);
  }
};

/**
 * Vectorized version of asin().
 * @param x Container of variables.
 * @tparam T Container type.
 * @return Arcsine of each variable in the container, in radians.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto asin(const T& x) {
  return apply_scalar_unary<asin_fun, T>::apply(x);
}

/**
 * Version of asin() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Arcsine of each variable in the container, in radians.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto asin(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().asin().matrix();
}

}  // namespace math
}  // namespace stan

#endif
