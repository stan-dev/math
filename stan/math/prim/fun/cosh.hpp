#ifndef STAN_MATH_PRIM_FUN_COSH_HPP
#define STAN_MATH_PRIM_FUN_COSH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap cosh() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic cosine of x.
 */
struct cosh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::cosh;
    return cosh(x);
  }
};

/**
 * Vectorized version of cosh().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Hyberbolic cosine of x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline typename apply_scalar_unary<cosh_fun, T>::return_t cosh(const T& x) {
  return apply_scalar_unary<cosh_fun, T>::apply(x);
}

/**
 * Version of cosh() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Hyberbolic cosine of x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto cosh(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().cosh().matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
