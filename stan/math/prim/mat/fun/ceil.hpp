#ifndef STAN_MATH_PRIM_MAT_FUN_CEIL_HPP
#define STAN_MATH_PRIM_MAT_FUN_CEIL_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap ceil() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Least integer >= x.
 */
struct ceil_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::ceil;
    return ceil(x);
  }
};

/**
 * Vectorized version of ceil().
 *
 * @tparam T type of container
 * @param x container
 * @return Least integer >= each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto ceil(const T& x) {
  return apply_scalar_unary<ceil_fun, T>::apply(x);
}

/**
 * Version of ceil() that accepts Eigen Matrix ar matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Least integer >= each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto ceil(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().ceil().matrix();
}

}  // namespace math
}  // namespace stan

#endif
