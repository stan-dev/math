#ifndef STAN_MATH_PRIM_PROB_GAMMA_POISSON_LPMF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_POISSON_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/neg_binomial_2_lpmf.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log-probability of a Poisson distribution with a gamma prior
 * for the rate parameter, after marginalizing the rate parameter out of the
 * likelihood. The likelihood is implemented using the neg_binomial_2_lpmf,
 * for more details on the derivation see:
 * https://gregorygundersen.com/blog/2019/09/16/poisson-gamma-nb/
 *
 *
 * @tparam propto
 * @tparam T_n Type of count outcome
 * @tparam T_shape Type of shape parameter for the Gamma prior
 * @tparam T_inv_scale Type of rate parameter for the Gamma prior
 * @param n Discrete count outcome
 * @param alpha Shape parameter for the Gamma prior
 * @param beta Rate parameter for the Gamma prior
 * @return log-likelihood
 */
template <bool propto, typename T_n, typename T_shape, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_shape, T_inv_scale>* = nullptr>
return_type_t<T_shape, T_inv_scale> gamma_poisson_lpmf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  static const char* function = "gamma_poisson_lpmf";
  // To avoid an integer division below, the shape parameter is promoted to a
  // double if it is an integer
  using AlphaScalarT = scalar_type_t<T_shape>;
  using PromotedIfIntT
      = std::conditional_t<std::is_integral<AlphaScalarT>::value, double,
                           AlphaScalarT>;

  const auto& n_ref = to_ref(n);
  ref_type_t<promote_scalar_t<PromotedIfIntT, T_shape>> alpha_ref = alpha;
  const auto& beta_ref = to_ref(beta);

  check_nonnegative(function, "Random variable", n_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  return neg_binomial_2_lpmf<propto>(n_ref, elt_divide(alpha_ref, beta_ref),
                                     alpha_ref);
}

template <typename T_n, typename T_shape, typename T_inv_scale>
inline return_type_t<T_shape, T_inv_scale> gamma_poisson_lpmf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  return gamma_poisson_lpmf<false>(n, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
