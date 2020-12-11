#ifndef STAN_MATH_PRIM_FUN_EXP_HPP
#define STAN_MATH_PRIM_FUN_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>
#include <limits>

namespace stan {
namespace math {
/**
 * Structure to wrap `exp()` so that it can be
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
 * Return the elementwise `exp()` of the specified argument,
 * which may be a scalar or any Stan container of numeric scalars.
 * The return type is the same as the argument type.
 *
 * @tparam Container type of container
 * @param[in] x container
 * @return Elementwise application of exponentiation to the argument.
 */
template <
    typename Container,
    require_not_container_st<std::is_arithmetic, Container>* = nullptr,
    require_not_nonscalar_prim_or_rev_kernel_expression_t<Container>* = nullptr,
    require_not_var_matrix_t<Container>* = nullptr>
inline auto exp(const Container& x) {
  return apply_scalar_unary<exp_fun, Container>::apply(x);
}

/**
 * Version of `exp()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Elementwise application of exponentiation to the argument.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto exp(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().exp(); });
}

namespace internal {
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
template <typename V>
inline std::complex<V> complex_exp(const std::complex<V>& z) {
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
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
