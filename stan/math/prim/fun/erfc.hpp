#ifndef STAN_MATH_PRIM_FUN_ERFC_HPP
#define STAN_MATH_PRIM_FUN_ERFC_HPP

#include <stan/math/prim/meta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the complementary error function of the specified value.
 *
 * \f[
 * \mbox{erfc}(x) = 1 - \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt
 * \f]
 *
 * @param[in] x Argument.
 * @return Complementary error function of the argument.
 */
inline double erfc(double x) { return std::erfc(x); }

/**
 * Return the error function of the specified argument.  This
 * version is required to disambiguate <code>erfc(int)</code>.
 *
 * @param[in] x Argument.
 * @return Complementary error function value of the argument.
 */
inline double erfc(int x) { return std::erfc(x); }

/**
 * Structure to wrap erfc() so that it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable
 * @return Complementary error function applied to x.
 */
struct erfc_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return erfc(x);
  }
};

/**
 * Vectorized version of erfc().
 *
 * @tparam T type of container
 * @param x container
 * @return Complementary error function applied to each value in x.
 */
template <typename T>
inline auto erfc(const T& x) {
  return apply_scalar_unary<erfc_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
