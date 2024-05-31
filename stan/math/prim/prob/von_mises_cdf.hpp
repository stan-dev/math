#ifndef STAN_MATH_PRIM_PROB_VON_MISES_CDF_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_CDF_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/modified_bessel_first_kind.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

/**
 * This implementation of the von Mises cdf
 * is a copy of scipy's. see:
 * https://github.com/scipy/scipy/blob/8dba340293fe20e62e173bdf2c10ae208286692f/scipy/stats/vonmises.py
 *
 * When k < 50, approximate the von Mises cdf
 * with the following expansion that comes from
 * scipy.
 */
template <typename T_x, typename T_k>
return_type_t<T_x, T_k> von_mises_cdf_series(const T_x& x, const T_k& k) {
  const double pi = stan::math::pi();
  int p = value_of_rec(28 + 0.5 * k - 100 / (k + 5) + 1);
  auto s = sin(x);
  auto c = cos(x);
  auto sn = sin(p * x);
  auto cn = cos(p * x);

  using return_t = return_type_t<T_x, T_k>;
  return_t R = 0.0;
  return_t V = 0.0;

  int n;
  for (n = 1; n < p; n++) {
    auto sn_tmp = sn * c - cn * s;
    cn = cn * c + sn * s;
    sn = sn_tmp;
    R = 1 / (2.0 * (p - n) / k + R);
    V = R * (sn / (p - n) + V);
  }
  return 0.5 + x / TWO_PI + V / pi;
}

/**
 * conv_mises_cdf_normapprox(x, k) is used to approximate the von
 * Mises cdf for k > 50. In this regime, the von Mises cdf
 * is well-approximated by a normal distribution.
 */
template <typename T_x, typename T_k>
return_type_t<T_x, T_k> von_mises_cdf_normalapprox(const T_x& x, const T_k& k) {
  using std::exp;
  using std::sqrt;

  const auto b = sqrt(2 / pi()) * exp(k) / modified_bessel_first_kind(0, k);
  const auto z = b * sin(x / 2);
  const double mu = 0;
  const double sigma = 1;

  return normal_cdf(z, mu, sigma);
}

/**
 * This function calculates the cdf of the von Mises distribution in
 * the case where the distribution has support on (-pi, pi) and
 * has mean 0. If k is sufficiently small, this function approximates
 * the cdf with a Gaussian. Otherwise, use the expansion from scipy.
 */
template <typename T_x, typename T_k>
return_type_t<T_x, T_k> von_mises_cdf_centered(const T_x& x, const T_k& k) {
  using return_t = return_type_t<T_x, T_k>;
  return_t f;
  if (k < 49) {
    f = von_mises_cdf_series(x, k);
    if (f < 0) {
      f = 0.0;
      return f;
    }
    if (f > 1) {
      f = 1.0;
      return f;
    }
    return f;
  } else if (k < 50) {
    f = (50.0 - k) * von_mises_cdf_series(x, 49.0)
        + (k - 49.0) * von_mises_cdf_normalapprox(x, 50.0);
    return f;
  } else {
    f = von_mises_cdf_normalapprox(x, k);
    return f;
  }
}

}  // namespace internal

/** \ingroup prob_dists
 * Calculates the cumulative distribution function of the von Mises
 * distribution:
 *
 * \f$VonMisesCDF(x, \mu, \kappa) = \frac{1}{2\pi I_0(\kappa)} \int_{-\pi}^x\f$
 * \f$e^{\kappa cos(t - \mu)} dt\f$
 *
 * where
 *
 * \f$x \in [-\pi, \pi]\f$, \f$\mu \in \mathbb{R}\f$,
 * and \f$\kappa \in \mathbb{R}^+\f$.
 *
 * @param x Variate on the interval \f$[-pi, pi]\f$
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
  using internal::von_mises_cdf_centered;
  const double pi = stan::math::pi();

  using T_return = return_type_t<T_x, T_mu, T_k>;
  using T_x_ref = ref_type_t<T_x>;
  using T_mu_ref = ref_type_t<T_mu>;
  using T_k_ref = ref_type_t<T_k>;
  static char const* function = "von_mises_cdf";
  check_consistent_sizes(function, "Random variable", x, "Location parameter",
                         mu, "Scale parameter", k);

  T_x_ref x_ref = x;
  T_mu_ref mu_ref = mu;
  T_k_ref k_ref = k;

  check_bounded(function, "Random variable", value_of(x_ref), -pi, pi);
  check_finite(function, "Location parameter", value_of(mu_ref));
  check_positive(function, "Scale parameter", value_of(k_ref));

  if (size_zero(x, mu, k)) {
    return 1.0;
  }

  T_return res(1.0);

  scalar_seq_view<T_x_ref> x_vec(x_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_k_ref> k_vec(k_ref);
  size_t N = max_size(x, mu, k);

  for (size_t n = 0; n < N; ++n) {
    auto x_n = x_vec[n];
    auto mu_n = mu_vec[n];
    auto k_n = k_vec[n];

    if (x_n == -pi) {
      res *= 0.0;
    } else if (x_n == pi) {
      res *= 1.0;
    } else {
      // shift x so that mean is 0
      T_return x2 = x_n - mu_n;

      // x2 is on an interval (2*n*pi, (2*n + 1)*pi), move it to (-pi, pi)
      x2 += pi;
      const auto x_floor = floor(x2 / TWO_PI);
      const auto x_modded = x2 - x_floor * TWO_PI;
      x2 = x_modded - pi;

      // mu is on an interval (2*n*pi, (2*n + 1)*pi), move it to (-pi, pi)
      T_return mu2 = mu_n + pi;
      const auto mu_floor = floor(mu2 / TWO_PI);
      const auto mu_modded = mu2 - mu_floor * TWO_PI;
      mu2 = mu_modded - pi;

      // find cut
      T_return cut, bound_val;
      if (mu2 < 0) {
        cut = mu2 + pi;
        bound_val = -pi - mu2;
      }
      if (mu2 >= 0) {
        cut = mu2 - pi;
        bound_val = pi - mu2;
      }

      T_return f_bound_val = von_mises_cdf_centered(bound_val, k_n);
      T_return cdf_n;
      if (x_n <= cut) {
        cdf_n = von_mises_cdf_centered(x2, k_n) - f_bound_val;
      } else {
        cdf_n = von_mises_cdf_centered(x2, k_n) + 1 - f_bound_val;
      }

      if (cdf_n < 0.0)
        cdf_n = 0.0;

      res *= cdf_n;
    }
  }

  return res;
}

}  // namespace math
}  // namespace stan

#endif
