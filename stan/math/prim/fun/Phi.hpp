#ifndef STAN_MATH_PRIM_FUN_PHI_HPP
#define STAN_MATH_PRIM_FUN_PHI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfcx.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {
namespace internal {

/**
 * erf(u) for |u| < 0.46875: the first Cody (1969) rational, no exp.
 *
 * @tparam T scalar type
 * @param u argument, |u| < 0.46875
 * @return erf(u)
 */
template <typename T>
inline T erf_cody_small(const T& u) {
  const T s = u * u;
  const T num = (((1.85777706184603153e-1 * s + 3.16112374387056560) * s
                  + 1.13864154151050156e2)
                     * s
                 + 3.77485237685302021e2)
                    * s
                + 3.20937758913846947e3;
  const T den = (((s + 2.36012909523441209e1) * s + 2.44024637934444173e2) * s
                 + 1.28261652607737228e3)
                    * s
                + 2.84423683343917062e3;
  return u * num / den;
}

}  // namespace internal

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
 * With `u = x / sqrt(2)`: `0.5 + 0.5 * erf(u)` for `|u| < 0.46875`, where
 * `Phi` is near 1/2; otherwise `0.5 * exp(-x^2 / 2) * erfcx(|u|)` for
 * `x < 0` and one minus it for `x > 0`. No branch cancels, and only `exp`
 * comes from the C library. `-(0.5 * x) * x` rounds once, so the lower tail
 * keeps its relative accuracy down to underflow.
 *
 * @param x Argument.
 * @return Probability random sample is less than or equal to argument.
 */
inline double Phi(double x) {
  check_not_nan("Phi", "x", x);
  const double u = INV_SQRT_TWO * x;
  const double a = std::fabs(u);
  if (a < 0.46875) {
    return 0.5 + 0.5 * internal::erf_cody_small(u);
  }
  const double tail = 0.5 * std::exp(-(0.5 * x) * x)
                      * (a >= 4.0 ? internal::erfcx_cody_tail(a)
                                  : internal::erfcx_cody_middle(a));
  return x < 0 ? tail : 1.0 - tail;
}

/**
 * Structure to wrap Phi() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument
 * @return Unit normal CDF of x.
 */
struct Phi_fun {
  template <typename T>
  static inline auto fun(T&& x) {
    return Phi(std::forward<T>(x));
  }
};

/**
 * Vectorized version of Phi().
 *
 * @tparam T type of container
 * @param x container
 * @return Unit normal CDF of each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_container_t<T>* = nullptr, require_not_var_matrix_t<T>* = nullptr>
inline auto Phi(T&& x) {
  return apply_scalar_unary<Phi_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
