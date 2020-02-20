#ifndef STAN_MATH_PRIM_FUN_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_SQUARE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the square of the specified argument.
 *
 * <p>\f$\mbox{square}(x) = x^2\f$.
 *
 * <p>The implementation of <code>square(x)</code> is just
 * <code>x * x</code>.  Given this, this method is mainly useful
 * in cases where <code>x</code> is not a simple primitive type,
 * particularly when it is an autodiff type.
 *
 * @param x Input to square.
 * @return Square of input.
 */
inline double square(double x) { return std::pow(x, 2); }

/**
 * Structure to wrap square() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return x squared.
 */
struct square_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return square(x);
  }
};

/**
 * Vectorized version of square().
 *
 * @tparam T type of container
 * @param x container
 * @return Each value in x squared.
 */
template <typename T,
          require_not_container_st<is_container, std::is_arithmetic, T>...>
inline auto square(const T& x) {
  return apply_scalar_unary<square_fun, T>::apply(x);
}

/**
 * Version of square() that accepts Eigen Matrix/Array objects or expressions.
 *
 * @tparam T Type of x
 * @param x Eigen Matrix/Array or expression
 * @return Each value in x squared.
 */
template <typename T,
          require_container_st<is_container, std::is_arithmetic, T>...>
inline auto square(const T& x) {
  return apply_vector_unary<T>::apply(
      x, [&](const auto& v) { return v.derived().array().square(); });
}

}  // namespace math
}  // namespace stan

#endif
