#ifndef STAN_MATH_PRIM_FUN_ACOS_HPP
#define STAN_MATH_PRIM_FUN_ACOS_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

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
 * @tparam Container type of container
 * @param x container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Container,
          require_not_container_st<is_container,
                                   std::is_arithmetic, Container>...>
inline auto acos(const Container& x) {
  return apply_scalar_unary<acos_fun, Container>::apply(x);
}

/**
 * Version of acos() that accepts std::vectors, Eigen Matrix/Array objects
 *  or expressions, and containers of these.
 *
 * @tparam Container Type of x
 * @param x Container
 * @return Arc cosine of each variable in the container, in radians.
 */
template <typename Container,
          require_container_st<is_container, std::is_arithmetic, Container>...>
inline auto acos(const Container& x) {
  return apply_vector_unary<Container>::apply(
      x, [](auto&& v) { return v.array().acos(); });
}

}  // namespace math
}  // namespace stan

#endif
