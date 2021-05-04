#ifndef STAN_MATH_PRIM_FUN_INV_PHI_HPP
#define STAN_MATH_PRIM_FUN_INV_PHI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The inverse of the unit normal cumulative distribution function.
 *
 * The return value for a specified input probability, $p$, is
 * the unit normal variate, $x$, such that
 *
 * \f$\Phi(x) = \int_{-\infty}^x \mbox{\sf Norm}(x|0, 1) \ dx = p\f$
 *
 * Algorithm first derived in 2003 by Peter Jon Aklam at
 * http://home.online.no/~pjacklam/notes/invnorm/
 *
 * @param p Argument between 0 and 1.
 * @return Real number
 */
inline double inv_Phi(double p) {
  check_bounded("inv_Phi", "Probability variable", p, 0, 1);

  if (p < 8e-311) {
    return NEGATIVE_INFTY;
  }
  if (p == 1) {
    return INFTY;
  }

  static const double a[6]
      = {-3.969683028665376e+01, 2.209460984245205e+02,  -2.759285104469687e+02,
         1.383577518672690e+02,  -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[5]
      = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[6]
      = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00,  2.938163982698783e+00};
  static const double d[4] = {7.784695709041462e-03, 3.224671290700398e-01,
                              2.445134137142996e+00, 3.754408661907416e+00};

  static const double p_low = 0.02425;
  static const double p_high = 0.97575;

  double x;
  if ((p_low <= p) && (p <= p_high)) {
    double q = p - 0.5;
    double r = square(q);
    x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5])
        * q
        / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
  } else if (p < p_low) {
    double q = std::sqrt(-2.0 * std::log(p));
    x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
        / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  } else {
    double q = std::sqrt(-2.0 * log1m(p));
    x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
        / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }

  if (x < 37.6) {  // gradient blows up past here
    double e = Phi(x) - p;
    double u = e * SQRT_TWO_PI * std::exp(0.5 * square(x));
    x -= u / (1.0 + 0.5 * x * u);
  }

  return x;
}

/**
 * Structure to wrap inv_Phi() so it can be vectorized.
 *
 * @tparam T type of variable
 * @param x variable in range [0, 1]
 * @return Inverse unit normal CDF of x.
 * @throw std::domain_error if x is not between 0 and 1.
 */
struct inv_Phi_fun {
  template <typename T>
  static inline T fun(const T& x) {
    return inv_Phi(x);
  }
};

/**
 * Vectorized version of inv_Phi().
 *
 * @tparam T type of container
 * @param x variables in range [0, 1]
 * @return Inverse unit normal CDF of each value in x.
 * @throw std::domain_error if any value is not between 0 and 1.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto inv_Phi(const T& x) {
  return apply_scalar_unary<inv_Phi_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
