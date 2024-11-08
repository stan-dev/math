#ifndef STAN_MATH_PRIM_FUN_TANH_HPP
#define STAN_MATH_PRIM_FUN_TANH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/core/operator_division.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the hyperbolic tangent of the arithmetic argument.
 *
 * @tparam V `Arithmetic` argument
 * @param[in] x argument
 * @return hyperbolic tangent of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto tanh(const T x) {
  return std::tanh(x);
}

/**
 * Return the hyperbolic tangent of the complex argument.
 *
 * @tparam V `complex<Arithmetic>` argument
 * @param[in] x argument
 * @return hyperbolic tangent of the argument
 */
template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto tanh(const T x) {
  return std::tanh(x);
}

/**
 * Structure to wrap `tanh()` so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic tangent of x.
 */
struct tanh_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return tanh(x);
  }
};

/**
 * Vectorized version of `tanh()`.
 *
 * @tparam Container type of container
 * @param x angles in radians
 * @return Hyperbolic tangent of each value in x.
 */
template <typename Container, require_ad_container_t<Container>* = nullptr>
inline auto tanh(const Container& x) {
  return apply_scalar_unary<tanh_fun, Container>::apply(x);
}

/**
 * Version of `tanh()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyperbolic tangent of each value in x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto tanh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().tanh(); });
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
  auto exp_z = exp(z);
  auto exp_neg_z = exp(-z);
  return stan::math::internal::complex_divide(exp_z - exp_neg_z,
                                              exp_z + exp_neg_z);
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
