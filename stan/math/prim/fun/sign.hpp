#ifndef STAN_MATH_PRIM_FUN_SIGN_HPP
#define STAN_MATH_PRIM_FUN_SIGN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

// returns 1 if NaN is passed in.
template <typename T, require_stan_scalar_t<T>* = nullptr>
inline int sign(const T& z) {
  return (z == 0) ? 0 : z < 0 ? -1 : 1;
}

/**
 * Structure to wrap `sign()` so it can be vectorized.
 */
struct sign_fun {
  /**
   * Return the sign of the specified argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return sign of the argument.
   */
  template <typename T>
  static inline int fun(const T& x) {
    return sign(x);
  }
};

/**
 * Return the elementwise application of `sign()` to
 * specified argument container.
 *
 * @tparam T type of container
 * @param x container
 * @return Elementwise sign of members of container.
 */
template <typename T, require_container_t<T>* = nullptr>
inline auto sign(const T& x) {
  return apply_scalar_unary<sign_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
