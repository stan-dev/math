#ifndef STAN_MATH_PRIM_PROB_VON_MISES_CDF_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_CDF_HPP

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/prob/normal_cdf.hpp>
#include <stan/math/prim/fun/modified_bessel_first_kind.hpp>

namespace stan {
namespace math {

namespace internal {

template <typename T_x, typename T_k, typename T_p>
return_type_t<T_x, T_k, T_p> von_mises_cdf_series(const T_x& x, const T_k& k,
                                                  const T_p& p) {
  using std::cos;
  using std::sin;

  double pi = stan::math::pi();
  auto s = sin(x);
  auto c = cos(x);
  auto sn = sin(p * x);
  auto cn = cos(p * x);

  double R = 0;
  return_type_t<T_x, T_k, T_p> V = 0;

  int n;
  for (n = p - 1; n > 0; n--) {
    auto sn_tmp = sn * c - cn * s;
    cn = cn * c + sn * s;
    sn = sn_tmp;
    R = 1 / (2.0 * n / k + R);
    V = R * (sn / n + V);
  }
  return 0.5 + x / (2 * pi) + V / pi;
}

template <typename T_x, typename T_k>
return_type_t<T_x, T_k> von_mises_cdf_normalapprox(const T_x& x, const T_k& k) {
  using std::exp;
  using std::sqrt;
  double pi = stan::math::pi();

  auto b = sqrt(2 / pi) * exp(k) / modified_bessel_first_kind(0, k);
  auto z = b * sin(x / 2);
  double mu = 0;
  double sigma = 1;

  return normal_cdf(z, mu, sigma);
}

template <typename T_x, typename T_k>
return_type_t<T_x, T_k> von_mises_cdf_centered(const T_x& x, const T_k& k) {
  // if the scale is sufficiently small, approximate the cdf with a
  // normal distribution, otherwise, use the expansion from scipy
  double ck = 50;
  double a[4] = {28, 0.5, 100, 5};
  return_type_t<T_x, T_k> f;
  if (k < ck) {
    int p = a[0] + a[1] * k - a[2] / (k + a[3]) + 1;
    f = von_mises_cdf_series(x, k, p);
    if (f < 0)
      f = 0;
    if (f > 1)
      f = 1;
  } else {
    f = von_mises_cdf_normalapprox(x, k);
  }
  return f;
}

}  // namespace internal
using namespace internal;

/** \ingroup prob_dists
 * Calculates the cumulative distribution function of the von Mises
 * distribution:
 *
 * \f$VonMisesCDF(x, mu, \kappa) = \frac{1}{2\pi I_0(\kappa)} \int_{-\mu-\pi}^x
 * e^{\kappa cos(t - \mu)} dt\f$
 *
 * where
 *
 * \f$\mu \in \mathbb{R}\f$, \f$x \in (\mu - \pi, \mu + \pi)\f$, and \f$\kappa
 * \in \mathbb{R}^+\f$.
 *
 * @param x A scalar variate
 * @param mu The mean of the distribution
 * @param k The inverse scale of the distriubtion
 * @return The von Mises cdf evaluated at the specified arguments
 * @tparam T_x Type of x
 * @tparam T_mu Type of mean parameter
 * @tparam T_k Type of inverse scale parameter
 */
template <typename T_x, typename T_mu, typename T_k>
inline return_type_t<T_x, T_mu, T_k> von_mises_cdf(const T_x& x, const T_mu& mu,
                                                   const T_k& k) {
  static char const* const function = "von_mises_cdf";
  using T_partials_return = partials_return_t<T_x, T_mu, T_k>;

  check_not_nan(function, "Random variable", x);
  check_not_nan(function, "Scale parameter", k);
  check_not_nan(function, "Location parameter", mu);
  check_positive(function, "Scale parameter", k);

  // shift x so that mean is 0
  auto x2 = x - mu;

  // x is on an interval (2*n*pi, (2*n + 1)*pi), move it to (-pi, pi)
  double pi = stan::math::pi();
  x2 += pi;
  auto x_floor = floor(x2 / (2 * pi));
  auto x_moded = x2 - x_floor * 2 * pi;
  x2 = x_moded - pi;

  return von_mises_cdf_centered(x2, k);
}

}  // namespace math
}  // namespace stan

#endif
