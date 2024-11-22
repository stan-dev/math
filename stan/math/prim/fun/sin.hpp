#ifndef STAN_MATH_PRIM_FUN_SIN_HPP
#define STAN_MATH_PRIM_FUN_SIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/i_times.hpp>
#include <stan/math/prim/fun/sinh.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

/**
 * Return the sine of the complex argument.
 *
 * @tparam V `Arithmetic` argument
 * @param[in] x argument
 * @return sine of the argument
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline auto sin(const T x) {
  return std::sin(x);
}

/**
 * Return the sine of the complex argument.
 *
 * @tparam V `complex<Arithmetic>` argument
 * @param[in] x argument
 * @return sine of the argument
 */
template <typename T, require_complex_bt<std::is_arithmetic, T>* = nullptr>
inline auto sin(const T x) {
  return std::sin(x);
}

/**
 * Structure to wrap sin() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Sine of x.
 */
struct sin_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return sin(x);
  }
};

/**
 * Vectorized version of sin().
 *
 * @tparam Container type of container
 * @param x angles in radians
 * @return Sine of each value in x.
 */
template <typename T, require_ad_container_t<T>* = nullptr>
inline auto sin(const T& x) {
  return apply_scalar_unary<sin_fun, T>::apply(x);
}

/**
 * Version of sin() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Sine of each value in x.
 */
template <typename Container,
          require_container_bt<std::is_arithmetic, Container>* = nullptr>
inline auto sin(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().sin(); });
}

namespace internal {
/**
 * Return the sine of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return sine of the argument
 */
template <typename V>
inline std::complex<V> complex_sin(const std::complex<V>& z) {
  using namespace stan::math;
  return neg_i_times(sinh(i_times(z)));
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
