#ifndef STAN_MATH_PRIM_FUN_ABS_HPP
#define STAN_MATH_PRIM_FUN_ABS_HPP

#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/hypot.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>
#include <complex>

namespace stan {
namespace math {

template <typename T, require_arithmetic_t<T>* = nullptr>
auto abs(T x) {
  return std::abs(x);
}

template <typename T, require_complex_t<T>* = nullptr>
auto abs(T x) {
  return hypot(x.real(), x.imag());
}

/**
 * Structure to wrap `abs()` so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x argument
 * @return Absolute value of variable.
 */
struct abs_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return fabs(x);
  }
};

/**
 * Returns the elementwise `abs()` of the input,
 * which may be a scalar or any Stan container of numeric scalars.
 *
 * @tparam Container type of container
 * @param x argument
 * @return Absolute value of each variable in the container.
 */
template <typename Container,
          require_not_container_st<std::is_arithmetic, Container>* = nullptr,
          require_not_var_matrix_t<Container>* = nullptr,
          require_not_stan_scalar_t<Container>* = nullptr>
inline auto abs(const Container& x) {
  return apply_scalar_unary<abs_fun, Container>::apply(x);
}

/**
 * Version of `abs()` that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x argument
 * @return Absolute value of each variable in the container.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto abs(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [&](const auto& v) { return v.array().abs(); });
}

namespace internal {
/**
 * Return the absolute value of the complex argument.
 *
 * @tparam V value type of argument
 * @param[in] z argument
 * @return absolute value of the argument
 */
template <typename V>
inline V complex_abs(const std::complex<V>& z) {
  return hypot(z.real(), z.imag());
}
}  // namespace internal

}  // namespace math
}  // namespace stan

#endif
