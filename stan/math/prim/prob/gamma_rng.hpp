#ifndef STAN_MATH_PRIM_PROB_GAMMA_RNG_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a gamma random variate for the given shape and inverse
 * scale parameters using the specified random number generator.
 *
 * alpha and beta can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_shape type of shape parameter
 * @tparam T_inv type of inverse scale parameter
 * @tparam RNG type of random number generator
 * @param alpha (Sequence of) positive shape parameter(s)
 * @param beta (Sequence of) positive inverse scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) gamma random variate(s)
 * @throw std::domain_error if alpha or beta are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_shape, typename T_inv, class RNG>
inline typename VectorBuilder<true, double, T_shape, T_inv>::type gamma_rng(
    const T_shape& alpha, const T_inv& beta, RNG& rng) {
  using boost::gamma_distribution;
  using boost::variate_generator;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv>;
  static constexpr const char* function = "gamma_rng";
  check_consistent_sizes(function, "Shape parameter", alpha,
                         "Inverse scale Parameter", beta);
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(alpha, beta);
  VectorBuilder<true, double, T_shape, T_inv> output(N);

  for (size_t n = 0; n < N; ++n) {
    // Convert rate (inverse scale) argument to scale for boost
    variate_generator<RNG&, gamma_distribution<> > gamma_rng(
        rng, gamma_distribution<>(alpha_vec[n],
                                  1 / static_cast<double>(beta_vec[n])));
    output[n] = gamma_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
