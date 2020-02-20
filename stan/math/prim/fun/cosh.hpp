#ifndef STAN_MATH_PRIM_FUN_COSH_HPP
#define STAN_MATH_PRIM_FUN_COSH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap cosh() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic cosine of x.
 */
struct cosh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::cosh;
    return cosh(x);
  }
};

/**
 * Vectorized version of cosh().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Hyberbolic cosine of x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto cosh(const T& x) {
  return apply_scalar_unary<cosh_fun, T>::apply(x);
}

/**
 * Version of cosh() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Hyberbolic cosine of x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto cosh(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().cosh();
  });
}
}  // namespace math
}  // namespace stan

#endif
