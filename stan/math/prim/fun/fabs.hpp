#ifndef STAN_MATH_PRIM_FUN_FABS_HPP
#define STAN_MATH_PRIM_FUN_FABS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap fabs() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Absolute value of x.
 */
struct fabs_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::fabs;
    return fabs(x);
  }
};

/**
 * Vectorized version of fabs().
 *
 * @tparam T type of container
 * @param x container
 * @return Absolute value of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline typename apply_scalar_unary<fabs_fun, T>::return_t fabs(const T& x) {
  return apply_scalar_unary<fabs_fun, T>::apply(x);
}

/**
 * Version of fabs() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Absolute value of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto fabs(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().abs().matrix().eval();
}

/**
 * Version of fabs() that accepts Eigen Array or array expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Absolute value of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto fabs(const Eigen::ArrayBase<Derived>& x) {
  return x.derived().abs().eval();
}

}  // namespace math
}  // namespace stan

#endif
