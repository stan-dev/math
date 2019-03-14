#ifndef STAN_MATH_PRIM_FUN_ERFC_HPP
#define STAN_MATH_PRIM_FUN_ERFC_HPP

#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/fun/erfc.hpp>



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
 * version is required to disambiguate <code>erf(int)</code>.
 *
 * @param[in] x Argument.
 * @return Complementary error function value of the argument.
 */
inline double erfc(int x) { return std::erfc(x); }

}  // namespace math
}  // namespace stan







namespace stan {
namespace math {

/**
 * Structure to wrap erfc() so that it can be vectorized.
 * @param x Variable.
 * @tparam T Variable type.
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
 * @param x Container.
 * @tparam T Container type.
 * @return Complementary error function applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<erfc_fun, T>::return_t erfc(const T& x) {
  return apply_scalar_unary<erfc_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
