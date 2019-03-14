#ifndef STAN_MATH_PRIM_FUN_SQUARE_HPP
#define STAN_MATH_PRIM_FUN_SQUARE_HPP

#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>

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
 * particularly when it is an auto-dif type.
 *
 * @param x Input to square.
 * @return Square of input.
 */
inline double square(double x) { return x * x; }

/**
 * Structure to wrap square() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
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
 * @param x Container.
 * @tparam T Container type.
 * @return Each value in x squared.
 */
template <typename T>
inline typename apply_scalar_unary<square_fun, T>::return_t square(const T& x) {
  return apply_scalar_unary<square_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
