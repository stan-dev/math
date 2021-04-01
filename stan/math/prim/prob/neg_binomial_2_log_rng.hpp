#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_RNG_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a negative binomial random variate with the specified log-location
 * and inverse dispersion parameters using the given random number generator.
 *
 * eta and phi can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_loc Type of log-location parameter
 * @tparam T_inv Type of inverse overdispersion parameter
 * @tparam RNG type of random number generator
 * @param eta (Sequence of) positive log-location parameter(s)
 * @param phi (Sequence of) positive inverse dispersion parameter(s)
 * @param rng random number generator
 * @return (Sequence of) negative binomial random variate(s)
 * @throw std::domain_error if eta is non-finite or phi is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_inv, class RNG>
inline typename VectorBuilder<true, int, T_loc, T_inv>::type
neg_binomial_2_log_rng(const T_loc& eta, const T_inv& phi, RNG& rng) {
  using boost::gamma_distribution;
  using boost::variate_generator;
  using boost::random::poisson_distribution;
  using T_eta_ref = ref_type_t<T_loc>;
  using T_phi_ref = ref_type_t<T_inv>;
  static const char* function = "neg_binomial_2_log_rng";
  check_consistent_sizes(function, "Log-location parameter", eta,
                         "Inverse dispersion parameter", phi);
  T_eta_ref eta_ref = eta;
  T_phi_ref phi_ref = phi;
  check_finite(function, "Log-location parameter", eta_ref);
  check_positive_finite(function, "Inverse dispersion parameter", phi_ref);

  scalar_seq_view<T_eta_ref> eta_vec(eta_ref);
  scalar_seq_view<T_phi_ref> phi_vec(phi_ref);
  size_t N = max_size(eta, phi);
  VectorBuilder<true, int, T_loc, T_inv> output(N);

  for (size_t n = 0; n < N; ++n) {
    double exp_eta_div_phi
        = std::exp(static_cast<double>(eta_vec[n])) / phi_vec[n];

    // gamma_rng params must be positive and finite
    check_positive_finite(function,
                          "Exponential of the log-location parameter "
                          "divided by the precision parameter",
                          exp_eta_div_phi);

    double rng_from_gamma = variate_generator<RNG&, gamma_distribution<> >(
        rng, gamma_distribution<>(phi_vec[n], exp_eta_div_phi))();

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
