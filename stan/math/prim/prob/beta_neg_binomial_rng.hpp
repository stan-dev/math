#ifndef STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_BETA_NEG_BINOMIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/neg_binomial_rng.hpp>
#include <stan/math/prim/prob/beta_rng.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a beta-negative binomial random variate with given number of
 * successes, prior success, and prior failure parameters, using the given
 * random number generator.
 *
 * r, alpha, and beta can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_r type of number of successes parameter
 * @tparam T_alpha type of prior success parameter
 * @tparam T_beta type of prior failure parameter
 * @tparam RNG type of random number generator
 *
 * @param r (Sequence of) number of successes parameter(s)
 * @param alpha (Sequence of) positive success parameter(s)
 * @param beta (Sequence of) positive failure parameter(s)
 * @param rng random number generator
 * @return (Sequence of) beta-binomial random variate(s)
 * @throw std::domain_error if r is negative, or alpha or beta are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different sizes
 */
template <typename T_r, typename T_alpha, typename T_beta, class RNG>
inline typename VectorBuilder<true, int, T_r, T_alpha, T_beta>::type
beta_neg_binomial_rng(const T_r &r, const T_alpha &alpha, const T_beta &beta,
                      RNG &rng) {
  using T_r_ref = ref_type_t<T_r>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_beta_ref = ref_type_t<T_beta>;
  static constexpr const char *function = "beta_neg_binomial_rng";
  check_consistent_sizes(function, "Number of successes parameter", r,
                         "Prior success parameter", alpha,
                         "Prior failure parameter", beta);

  T_r_ref r_ref = r;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Number of successes parameter", r_ref);
  check_positive_finite(function, "Prior success parameter", alpha_ref);
  check_positive_finite(function, "Prior failure parameter", beta_ref);

  using T_p = decltype(beta_rng(alpha_ref, beta_ref, rng));
  T_p p = beta_rng(alpha_ref, beta_ref, rng);

  scalar_seq_view<T_p> p_vec(p);
  size_t size_p = stan::math::size(p);
  VectorBuilder<true, double, T_p> odds_ratio_p(size_p);
  for (size_t n = 0; n < size_p; ++n) {
    odds_ratio_p[n] = stan::math::exp(stan::math::log(p_vec.val(n))
                                      - stan::math::log((1 - p_vec.val(n))));
  }

  return neg_binomial_rng(r_ref, odds_ratio_p.data(), rng);
}

}  // namespace math
}  // namespace stan
#endif
