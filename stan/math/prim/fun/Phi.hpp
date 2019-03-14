#ifndef STAN_MATH_PRIM_FUN_PHI_HPP
#define STAN_MATH_PRIM_FUN_PHI_HPP

#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/vectorize/apply_scalar_unary.hpp>
#include <stan/math/prim/fun/Phi.hpp>









namespace stan {
namespace math {

/**
 * The unit normal cumulative distribution function.
 *
 * The return value for a specified input is the probability that
 * a random unit normal variate is less than or equal to the
 * specified value, defined by
 *
 * \f$\Phi(x) = \int_{-\infty}^x \mbox{\sf Norm}(x|0, 1) \ dx\f$
 *
 * This function can be used to implement the inverse link function
 * for probit regression.
 *
 * Phi will underflow to 0 below -37.5 and overflow to 1 above 8
 *
 * @param x Argument.
 * @return Probability random sample is less than or equal to argument.
 */
inline double Phi(double x) {
  check_not_nan("Phi", "x", x);
  if (x < -37.5)
    return 0;
  else if (x < -5.0)
    return 0.5 * erfc(-INV_SQRT_2 * x);
  else if (x > 8.25)
    return 1;
  else
    return 0.5 * (1.0 + erf(INV_SQRT_2 * x));
}

}  // namespace math
}  // namespace stan







namespace stan {
namespace math {

/**
 * Structure to wrap Phi() so it can be vectorized.
 * @param x Argument variable.
 * @tparam T Argument type.
 * @return Unit normal CDF of x.
 */
struct Phi_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return Phi(x);
  }
};

/**
 * Vectorized version of Phi().
 * @param x Container.
 * @tparam T Container type.
 * @return Unit normal CDF of each value in x.
 */
template <typename T>
inline typename apply_scalar_unary<Phi_fun, T>::return_t Phi(const T& x) {
  return apply_scalar_unary<Phi_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
