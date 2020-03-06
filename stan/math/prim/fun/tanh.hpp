#ifndef STAN_MATH_PRIM_FUN_TANH_HPP
#define STAN_MATH_PRIM_FUN_TANH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/core/operator_division.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Structure to wrap tanh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic tangent of x.
 */
struct tanh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::tanh;
    return tanh(x);
  }
};

/**
 * Vectorized version of tanh().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Hyperbolic tangent of each value in x.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto tanh(const T& x) {
  return apply_scalar_unary<tanh_fun, T>::apply(x);
}

/**
 * Version of tanh() that accepts Eigen Matrix or matrix expressions.
 *
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Hyperbolic tangent of each value in x.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto tanh(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().tanh().matrix().eval();
}

namespace internal {
/**
 * Return the hyperbolic tangent of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return hyperbolic tangent of the argument
 */
template <typename V>
inline std::complex<V> complex_tanh(const std::complex<V>& z) {
  using std::exp;
  using stan::math::operator/;
  auto exp_z = exp(z);
  auto exp_neg_z = exp(-z);
  return stan::math::internal::complex_divide(exp_z - exp_neg_z, exp_z + exp_neg_z);
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
