#ifndef STAN_MATH_PRIM_FUN_ACOS_HPP
#define STAN_MATH_PRIM_FUN_ACOS_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/asin.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap acos() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Arc cosine of variable in radians.
 */
struct acos_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::acos;
    return acos(x);
  }
};

/**
 * Vectorized version of acos().
 *
 * @tparam T type of container
 * @param x container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto acos(const T& x) {
  return apply_scalar_unary<acos_fun, T>::apply(x);
}

/**
 * Version of acos() that accepts Eigen Matrix or matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto acos(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().acos().matrix().eval();
}

namespace internal {
/**
 * Return the arc cosine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return arc cosine of the argument
 */
template <typename V>
inline std::complex<V> complex_acos(const std::complex<V>& z) {
  return V(0.5 * pi()) - asin(z);
}
}



}  // namespace math
}  // namespace stan

#endif
