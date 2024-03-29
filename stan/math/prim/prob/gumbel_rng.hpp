#ifndef STAN_MATH_PRIM_PROB_GUMBEL_RNG_HPP
#define STAN_MATH_PRIM_PROB_GUMBEL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Gumbel random variate with the given location and scale
 * using the specified random number generator.
 *
 * mu and beta can each be a scalar or a vector. Any non-scalar inputs
 * must be the same length.
 *
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @tparam RNG type of random number generator
 * @param mu (Sequence of) location parameter(s)
 * @param beta (Sequence of) scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Gumbel random variate(s)
 * @throw std::domain_error if mu is infinite or beta is nonpositive.
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_scale, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_scale>::type gumbel_rng(
    const T_loc& mu, const T_scale& beta, RNG& rng) {
  using boost::uniform_01;
  using boost::variate_generator;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_beta_ref = ref_type_t<T_scale>;
  static constexpr const char* function = "gumbel_rng";
  check_consistent_sizes(function, "Location parameter", mu, "Scale Parameter",
                         beta);
  T_mu_ref mu_ref = mu;
  T_beta_ref beta_ref = beta;
  check_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Scale parameter", beta_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(mu, beta);
  VectorBuilder<true, double, T_loc, T_scale> output(N);

  variate_generator<RNG&, uniform_01<> > uniform01_rng(rng, uniform_01<>());
  for (size_t n = 0; n < N; ++n) {
    output[n] = mu_vec[n] - beta_vec[n] * std::log(-std::log(uniform01_rng()));
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
