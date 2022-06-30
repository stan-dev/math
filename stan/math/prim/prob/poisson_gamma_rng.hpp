#ifndef STAN_MATH_PRIM_PROB_POISSON_GAMMA_RNG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_GAMMA_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/neg_binomial_2_rng.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a poisson random variate where the rate parameter follows the
 * specified shape and inverse scale parameters using the given
 * random number generator.
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
 * @return (Sequence of) poisson-gamma random variate(s)
 * @throw std::domain_error if alpha or beta are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_shape, typename T_inv, class RNG>
inline typename VectorBuilder<true, int, T_shape, T_inv>::type
  poisson_gamma_rng(
    const T_shape& alpha, const T_inv& beta, RNG& rng) {
  static const char* function = "poisson_gamma_rng";
  // To avoid an integer division below, the shape parameter is promoted to a
  // double if it is an integer
  using AlphaScalarT = scalar_type_t<T_shape>;
  using PromotedIfIntT = std::conditional_t<
    std::is_integral<AlphaScalarT>::value,
    double, AlphaScalarT>;

  check_consistent_sizes(function, "Shape parameter", alpha,
                         "Inverse scale Parameter", beta);

  ref_type_t<promote_scalar_t<PromotedIfIntT, T_shape>> alpha_ref = alpha;
  const auto& beta_ref = to_ref(beta);

  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  return neg_binomial_2_rng(elt_divide(alpha_ref, beta_ref), alpha_ref, rng);
}

}  // namespace math
}  // namespace stan
#endif
