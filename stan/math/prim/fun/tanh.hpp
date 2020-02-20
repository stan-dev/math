#ifndef STAN_MATH_PRIM_FUN_TANH_HPP
#define STAN_MATH_PRIM_FUN_TANH_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Structure to wrap tanh() so that it can be vectorized.
 *
 * @tparam T type of argument
 * @param x angle in radians
 * @return Hyperbolic tangent of x.
 */
struct tanh_fun {
  template <typename T>
  static inline T fun(const T& x) {
    using std::tanh;
    return tanh(x);
  }
};

/**
 * Vectorized version of tanh().
 *
 * @tparam T type of container
 * @param x angles in radians
 * @return Hyperbolic tangent of each value in x.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto tanh(const T& x) {
  return apply_scalar_unary<tanh_fun, T>::apply(x);
}

/**
 * Version of tanh() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Hyperbolic tangent of each value in x.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto tanh(const T& x) {
  return apply_vector_unary<T>::apply(x, [&](const auto& v) {
    return v.derived().array().tanh();
  });
}

}  // namespace math
}  // namespace stan

#endif
