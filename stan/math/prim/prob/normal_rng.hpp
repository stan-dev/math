#ifndef STAN_MATH_PRIM_PROB_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Normal random variate for the given location and scale
 * using the specified random number generator.
 *
 * mu and sigma can each be a scalar or a vector. Any non-scalar
 * inputs must be the same length.
 *
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @tparam RNG type of random number generator
 * @param mu (Sequence of) location parameter(s)
 * @param sigma (Sequence of) positive scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) normal random variate(s)
 * @throw std::domain_error if mu is infinite or sigma is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_scale, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_scale>::type normal_rng(
    const T_loc& mu, const T_scale& sigma, RNG& rng) {
  using boost::normal_distribution;
  using boost::variate_generator;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static const char* function = "normal_rng";
  check_consistent_sizes(function, "Location parameter", mu, "Scale Parameter",
                         sigma);
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(mu, sigma);
  VectorBuilder<true, double, T_loc, T_scale> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, normal_distribution<> > norm_rng(
        rng, normal_distribution<>(mu_vec[n], sigma_vec[n]));
    output[n] = norm_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
