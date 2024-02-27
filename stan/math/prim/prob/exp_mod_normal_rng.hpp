#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/prob/exponential_rng.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/prob/normal_rng.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return an exponentially modified normal random variate for the
 * given location, scale, and inverse scale using the specified random
 * number generator.
 *
 * mu, sigma, and lambda can each be a scalar or a one-dimensional container.
 * Any non-scalar inputs must be the same size.
 *
 * @tparam T_loc Type of location parameter
 * @tparam T_scale Type of scale parameter
 * @tparam T_inv_scale Type of inverse scale parameter
 * @tparam RNG type of random number generator
 * @param mu (Sequence of) location parameter(s)
 * @param sigma (Sequence of) scale parameter(s)
 * @param lambda (Sequence of) inverse scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Exponentially modified normal random variate(s)
 * @throw std::domain_error if mu is infinite, sigma is nonpositive,
 * or lambda is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_scale, typename T_inv_scale, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_scale, T_inv_scale>::type
exp_mod_normal_rng(const T_loc& mu, const T_scale& sigma,
                   const T_inv_scale& lambda, RNG& rng) {
  static constexpr const char* function = "exp_mod_normal_rng";
  using T_mu_ref = ref_type_t<T_loc>;
  using T_sigma_ref = ref_type_t<T_scale>;
  using T_lambda_ref = ref_type_t<T_inv_scale>;
  check_consistent_sizes(function, "Location parameter", mu, "Scale Parameter",
                         sigma, "Inv_scale Parameter", lambda);

  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_lambda_ref lambda_ref = lambda;
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);
  check_positive_finite(function, "Inv_scale parameter", lambda_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  scalar_seq_view<T_lambda_ref> lambda_vec(lambda_ref);
  size_t N = max_size(mu, sigma, lambda);
  VectorBuilder<true, double, T_loc, T_scale, T_inv_scale> output(N);

  for (size_t n = 0; n < N; ++n) {
    output[n] = normal_rng(mu_vec[n], sigma_vec[n], rng)
                + exponential_rng(lambda_vec[n], rng);
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
