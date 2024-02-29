#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/normal_rng.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Skew-normal random variate for the given location, scale,
 * and shape using the specified random number generator.
 *
 * mu, sigma, and alpha can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @tparam T_shape type of shape parameter
 * @tparam RNG type of random number generator
 *
 * @param mu (Sequence of) location parameter(s)
 * @param sigma (Sequence of) scale parameter(s)
 * @param alpha (Sequence of) shape parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Skew-normal random variate(s)
 * @throw std::domain_error if mu is infinite, sigma is nonpositive, or
 * alpha is infinite
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_scale, typename T_shape, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_scale, T_shape>::type
skew_normal_rng(const T_loc& mu, const T_scale& sigma, const T_shape& alpha,
                RNG& rng) {
  using boost::variate_generator;
  using boost::random::normal_distribution;
  static constexpr const char* function = "skew_normal_rng";
  check_consistent_sizes(function, "Location parameter", mu, "Scale Parameter",
                         sigma, "Shape Parameter", alpha);
  const auto& mu_ref = to_ref(mu);
  const auto& sigma_ref = to_ref(sigma);
  const auto& alpha_ref = to_ref(alpha);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);
  check_finite(function, "Shape parameter", alpha_ref);

  scalar_seq_view<T_loc> mu_vec(mu_ref);
  scalar_seq_view<T_scale> sigma_vec(sigma_ref);
  scalar_seq_view<T_shape> alpha_vec(alpha_ref);
  size_t N = max_size(mu, sigma, alpha);
  VectorBuilder<true, double, T_loc, T_scale, T_shape> output(N);

  variate_generator<RNG&, normal_distribution<> > norm_rng(
      rng, normal_distribution<>(0, 1));
  for (size_t n = 0; n < N; ++n) {
    double r1 = norm_rng();
    double r2 = norm_rng();

    if (r2 > alpha_vec[n] * r1) {
      r1 = -r1;
    }

    output[n] = mu_vec[n] + sigma_vec[n] * r1;
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
