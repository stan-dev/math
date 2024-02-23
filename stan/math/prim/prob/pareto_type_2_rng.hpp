#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_RNG_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/prob/exponential_rng.hpp>
#include <stan/math/prim/prob/normal_rng.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Pareto type 2 random variate for the given location,
 * scale, and shape using the specified random number generator.
 *
 * mu, lambda, and alpha can each be a scalar or a one-dimensional container.
 * Any non-scalar inputs must be the same size.
 *
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @tparam T_shape type of shape parameter
 * @tparam RNG type of random number generator
 *
 * @param mu (Sequence of) location parameter(s)
 * @param lambda (Sequence of) scale parameter(s)
 * @param alpha (Sequence of) shape parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Pareto type 2 random variate(s)
 * @throw std::domain_error if mu is infinite or lambda or alpha are
 * nonpositive,
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_scale, typename T_shape, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_scale, T_shape>::type
pareto_type_2_rng(const T_loc& mu, const T_scale& lambda, const T_shape& alpha,
                  RNG& rng) {
  using boost::variate_generator;
  using boost::random::uniform_real_distribution;
  static constexpr const char* function = "pareto_type_2_rng";
  check_consistent_sizes(function, "Location parameter", mu, "Scale Parameter",
                         lambda, "Shape Parameter", alpha);
  const auto& mu_ref = to_ref(mu);
  const auto& lambda_ref = to_ref(lambda);
  const auto& alpha_ref = to_ref(alpha);
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", lambda_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);

  scalar_seq_view<T_loc> mu_vec(mu_ref);
  scalar_seq_view<T_scale> lambda_vec(lambda_ref);
  scalar_seq_view<T_shape> alpha_vec(alpha_ref);
  size_t N = max_size(mu, lambda, alpha);
  VectorBuilder<true, double, T_loc, T_scale, T_shape> output(N);

  variate_generator<RNG&, uniform_real_distribution<> > uniform_rng(
      rng, uniform_real_distribution<>(0.0, 1.0));
  for (size_t n = 0; n < N; ++n) {
    output[n] = (std::pow(1.0 - uniform_rng(), -1.0 / alpha_vec[n]) - 1.0)
                    * lambda_vec[n]
                + mu_vec[n];
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
