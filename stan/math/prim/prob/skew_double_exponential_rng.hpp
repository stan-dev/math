#ifndef STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a skew double exponential random variate with the given location
 * scale and skewness using the specified random number generator.
 *
 * mu, sigma and tau can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_loc Type of location parameter
 * @tparam T_scale Type of scale parameter
 * @tparam T_skewness Type of skewness parameter
 * @tparam RNG class of random number generator
 * @param mu (Sequence of) location parameter(s)
 * @param sigma (Sequence of) scale parameter(s)
 * @param tau (Sequence of) skewness parameter(s)
 * @param rng random number generator
 * @return (Sequence of) double exponential random variate(s)
 * @throw std::domain_error if mu is infinite or sigma is nonpositive or tau is
 *  not bound between 0.0 and 1.0
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_scale, typename T_skewness, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_scale, T_skewness>::type
skew_double_exponential_rng(const T_loc& mu, const T_scale& sigma,
                            const T_skewness& tau, RNG& rng) {
  using boost::variate_generator;
  using boost::random::uniform_real_distribution;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using T_tau_ref = ref_type_t<T_skewness>;
  static const char* function = "skew_double_exponential_rng";
  check_consistent_sizes(function, "Location parameter", mu, "Scale Parameter",
                         sigma, "Skewness Parameter", tau);
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_tau_ref tau_ref = tau;
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);
  check_bounded(function, "Skewness parameter", tau_ref, 0.0, 1.0);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  scalar_seq_view<T_tau_ref> tau_vec(tau_ref);
  size_t N = max_size(mu, sigma, tau);
  VectorBuilder<true, double, T_loc, T_scale, T_skewness> output(N);

  variate_generator<RNG&, uniform_real_distribution<> > z_rng(
      rng, uniform_real_distribution<>(0.0, 1.0));
  for (size_t n = 0; n < N; ++n) {
    double z = z_rng();
    if (z < tau_vec[n]) {
      output[n]
          = log(z / tau_vec[n]) * sigma_vec[n] / (2.0 * (1.0 - tau_vec[n]))
            + mu_vec[n];
    } else {
      output[n] = log((1.0 - z) / (1.0 - tau_vec[n])) * (-sigma_vec[n])
                      / (2.0 * tau_vec[n])
                  + mu_vec[n];
    }
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
