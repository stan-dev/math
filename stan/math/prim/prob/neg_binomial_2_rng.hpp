#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_RNG_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a negative binomial random variate with the specified location and
 * precision parameters using the given random number generator.
 *
 * mu and phi can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_loc type of location parameter
 * @tparam T_prec type of precision parameter
 * @tparam RNG type of random number generator
 * @param mu (Sequence of) positive location parameter(s)
 * @param phi (Sequence of) positive precision parameter(s)
 * @param rng random number generator
 * @return (Sequence of) negative binomial random variate(s)
 * @throw std::domain_error if mu or phi are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_prec, class RNG>
inline typename VectorBuilder<true, int, T_loc, T_prec>::type
neg_binomial_2_rng(const T_loc& mu, const T_prec& phi, RNG& rng) {
  using boost::gamma_distribution;
  using boost::variate_generator;
  using boost::random::poisson_distribution;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_phi_ref = ref_type_t<T_prec>;
  static constexpr const char* function = "neg_binomial_2_rng";
  check_consistent_sizes(function, "Location parameter", mu,
                         "Precision parameter", phi);
  T_mu_ref mu_ref = mu;
  T_phi_ref phi_ref = phi;
  check_positive_finite(function, "Location parameter", mu_ref);
  check_positive_finite(function, "Precision parameter", phi_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t N = max_size(mu, phi);
  VectorBuilder<true, int, T_loc, T_prec> output(N);

  for (size_t n = 0; n < N; ++n) {
    double mu_div_phi = static_cast<double>(mu_vec[n]) / phi_vec[n];

    // gamma_rng params must be positive and finite
    check_positive_finite(function,
                          "Location parameter divided by the "
                          "precision parameter",
                          mu_div_phi);

    double rng_from_gamma = variate_generator<RNG&, gamma_distribution<> >(
        rng, gamma_distribution<>(phi_vec[n], mu_div_phi))();

    // same as the constraints for poisson_rng
    check_less(function, "Random number that came from gamma distribution",
               rng_from_gamma, POISSON_MAX_RATE);
    check_not_nan(function, "Random number that came from gamma distribution",
                  rng_from_gamma);
    check_nonnegative(function,
                      "Random number that came from gamma distribution",
                      rng_from_gamma);

    output[n] = variate_generator<RNG&, poisson_distribution<> >(
        rng, poisson_distribution<>(rng_from_gamma))();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
