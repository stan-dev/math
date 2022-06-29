#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/prob/neg_binomial_2_rng.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a negative binomial random variate with the specified shape and
 * inverse scale parameters using the given random number generator.
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
 * @return (Sequence of) negative binomial random variate(s)
 * @throw std::domain_error if alpha or beta are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_shape, typename T_inv, class RNG>
inline typename VectorBuilder<true, int, T_shape, T_inv>::type poisson_gamma_rng(
    const T_shape& alpha, const T_inv& beta, RNG& rng) {
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv>;
  static const char* function = "poisson_gamma_rng";
  check_consistent_sizes(function, "Shape parameter", alpha,
                         "Inverse scale Parameter", beta);
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  return neg_binomial_2_rng(elt_divide(alpha_ref, beta_ref), alpha_ref);
}

}  // namespace math
}  // namespace stan
#endif
