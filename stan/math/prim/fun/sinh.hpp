#ifndef STAN_MATH_PRIM_FUN_SINH_HPP
#define STAN_MATH_PRIM_FUN_SINH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap sinh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic sine of x.
 */
struct sinh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::sinh;
    return sinh(x);
  }
};

/**
 * Vectorized version of sinh().
 *
 * @tparam Container type of container
 * @param x container
 * @return Hyperbolic sine of each variable in x.
 */
template <
    typename Container,
    require_not_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto sinh(const Container& x) {
  return apply_scalar_unary<sinh_fun, Container>::apply(x);
}

/**
 * Version of sinh() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Hyperbolic sine of each variable in x.
 */
template <typename Container,
          require_container_st<std::is_arithmetic, Container>* = nullptr>
inline auto sinh(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](const auto& v) { return v.array().sinh(); });
}

}  // namespace math
}  // namespace stan

#endif
