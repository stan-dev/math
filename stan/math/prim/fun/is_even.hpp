#ifndef STAN_MATH_PRIM_FUN_IS_EVEN_HPP
#define STAN_MATH_PRIM_FUN_IS_EVEN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * Return 1 if the integer is even and 0 otherwise.
 *
 * @tparam T Integer type
 * @param n Integer to check
 * @return 1 if n is even, 0 otherwise
 */
template <typename T, require_integral_t<T>* = nullptr>
inline int is_even(T n) {
  return n % 2 == 0 ? 1 : 0;
}

/**
 * Vectorization wrapper for is_even.
 */
struct is_even_fun {
  template <typename T>
  static inline int fun(const T& x) {
    return is_even(static_cast<int>(x));
  }
};

/**
 * Apply is_even elementwise to a container.
 *
 * @tparam T Container type
 * @param x Container of integers
 * @return Elementwise parity results
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto is_even(T&& x) {
  return apply_scalar_unary<is_even_fun, T>::apply(
      std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif