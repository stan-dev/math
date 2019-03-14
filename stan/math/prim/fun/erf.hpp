#ifndef STAN_MATH_PRIM_FUN_ERF_HPP
#define STAN_MATH_PRIM_FUN_ERF_HPP

#include <cmath>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/fun/erf.hpp>



namespace stan {
namespace math {

/**
 * Return the error function of the specified value.
 *
 * \f[
 * \mbox{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt
 * \f]
 *
 * @param[in] x Argument.
 * @return Error function of the argument.
 */
inline double erf(double x) { return std::erf(x); }

/**
 * Return the error function of the specified argument.  This
 * version is required to disambiguate <code>erf(int)</code>.
 *
 * @param[in] x Argument.
 * @return Error function of the argument.
 */
inline double erf(int x) { return std::erf(x); }

}  // namespace math
}  // namespace stan







namespace stan {
namespace math {

/**
 * Structure to wrap erf() so it can be vectorized.
 *
 * @tparam T Variable type.
 * @param x Variable.
 * @return Error function of x.
 */
struct erf_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return erf(x);
  }
};

/**
 * Vectorized version of erf().
 *
 * @tparam T Container type.
 * @param x Container.
 * @return Error function applied to each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<erf_fun, T>::return_t erf(const T& x) {
  return apply_scalar_unary<erf_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
