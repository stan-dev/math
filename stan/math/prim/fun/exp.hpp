#ifndef STAN_MATH_PRIM_FUN_EXP_HPP
#define STAN_MATH_PRIM_FUN_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the natural exponential of the specified argument.  This
 * version is required to disambiguate <code>exp(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural exponential of argument.
 */
inline double exp(int x) { return std::exp(x); }

/**
 * Return the natural (base e) complex exponentiation of the specified
 * complex argument.
 *
 * @tparam V value type (must be Stan autodiff type)
 * @param z complex number
 * @return natural exponentiation of specified complex number
 * @see documentation for `std::complex` for boundary condition and
 * branch cut details
 */
template <typename V, require_autodiff_t<V>>
inline std::complex<V> exp(const std::complex<V>& z) {
  if (is_inf(z.real()) && z.real() > 0) {
    if (is_nan(z.imag()) || z.imag() == 0) {
      // (+inf, nan), (+inf, 0)
      return z;
    } else if (is_inf(z.imag()) && z.imag() > 0) {
      // (+inf, +inf)
      return {z.real(), std::numeric_limits<double>::quiet_NaN()};
    } else if (is_inf(z.imag()) && z.imag() < 0) {
      // (+inf, -inf)
      return {std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};
    }
  }
  if (is_inf(z.real()) && z.real() < 0
      && (is_nan(z.imag()) || is_inf(z.imag()))) {
    // (-inf, nan), (-inf, -inf), (-inf, inf)
    return {0, 0};
  }
  if (is_nan(z.real()) && z.imag() == -0.0) {
    // (nan, -0)
    return z;
  }
  V exp_re = exp(z.real());
  return {exp_re * cos(z.imag()), exp_re * sin(z.imag())};
}

/**
 * Structure to wrap <code>exp()</code> so that it can be
 * vectorized.
 */
struct exp_fun {
  /**
   * Return the exponential of the specified scalar argument.
   *
   * @tparam T type of argument
   * @param[in] x argument
   * @return Exponential of argument.
   */
  template <typename T>
  static inline T fun(const T& x) {
    using std::exp;
    return exp(x);
  }
};

/**
 * Return the elementwise exponentiation of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam T type of container
 * @param[in] x container
 * @return Elementwise application of exponentiation to the argument.
 */
template <typename T, typename = require_not_eigen_vt<std::is_arithmetic, T>>
inline auto exp(const T& x) {
  return apply_scalar_unary<exp_fun, T>::apply(x);
}

/**
 * Version of exp() that accepts Eigen Matrix or matrix expressions.
 * @tparam Derived derived type of x
 * @param x Matrix or matrix expression
 * @return Elementwise application of exponentiation to the argument.
 */
template <typename Derived,
          typename = require_eigen_vt<std::is_arithmetic, Derived>>
inline auto exp(const Eigen::MatrixBase<Derived>& x) {
  return x.derived().array().exp().matrix().eval();
}

}  // namespace math
}  // namespace stan

#endif
