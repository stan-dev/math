#ifndef STAN_MATH_PRIM_PROB_VON_MISES_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_LCCDF_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/prob/von_mises_cdf.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the log of the complement of the cumulative distribution
 * function of the von Mises distribution:
 *
 * \f$VonMisesLCCDF(x, \mu, \kappa) = \log ( 1.0 - \frac{1}{2\pi I_0(\kappa)}
 * \int_{-\pi}^x\f$ \f$e^{\kappa cos(t - \mu)} dt )\f$
 *
 * where
 *
 * \f$x \in [-\pi, \pi]\f$, \f$\mu \in \mathbb{R}\f$,
 * and \f$\kappa \in \mathbb{R}^+\f$.
 *
 * @param x Variate on the interval \f$[-pi, pi]\f$
 * @param mu The mean of the distribution
 * @param k The inverse scale of the distriubtion
 * @return The log of the von Mises cdf evaluated at the specified arguments
 * @tparam T_x Type of x
 * @tparam T_mu Type of mean parameter
 * @tparam T_k Type of inverse scale parameter
 */
template <typename T_x, typename T_mu, typename T_k>
inline return_type_t<T_x, T_mu, T_k> von_mises_lccdf(const T_x& x,
                                                     const T_mu& mu,
                                                     const T_k& k) {
  using internal::von_mises_cdf_centered;
  const double pi = stan::math::pi();

  using T_return = return_type_t<T_x, T_mu, T_k>;
  using T_x_ref = ref_type_t<T_x>;
  using T_mu_ref = ref_type_t<T_mu>;
  using T_k_ref = ref_type_t<T_k>;
  static char const* function = "von_mises_lccdf";
  check_consistent_sizes(function, "Random variable", x, "Location parameter",
                         mu, "Scale parameter", k);

  T_x_ref x_ref = x;
  T_mu_ref mu_ref = mu;
  T_k_ref k_ref = k;

  check_bounded(function, "Random variable", value_of(x_ref), -pi, pi);
  check_finite(function, "Location parameter", value_of(mu_ref));
  check_positive(function, "Scale parameter", value_of(k_ref));

  if (size_zero(x, mu, k)) {
    return 0.0;
  }

  T_return res(0.0);

  scalar_seq_view<T_x_ref> x_vec(x_ref);
  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_k_ref> k_vec(k_ref);
  size_t N = max_size(x, mu, k);

  for (size_t n = 0; n < N; ++n) {
    auto x_n = x_vec[n];
    auto mu_n = mu_vec[n];
    auto k_n = k_vec[n];

    if (x_n == -pi) {
      res += 0.0;
    } else if (x_n == pi) {
      res += NEGATIVE_INFTY;
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

      res += log1m(cdf_n);
    }
  }

  return res;
}

}  // namespace math
}  // namespace stan

#endif
