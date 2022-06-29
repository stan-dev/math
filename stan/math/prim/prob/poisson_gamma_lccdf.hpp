#ifndef STAN_MATH_PRIM_PROB_POISSON_GAMMA_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_GAMMA_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/prob/neg_binomial_2_lccdf.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_n, typename T_shape, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_shape, T_inv_scale>* = nullptr>
return_type_t<T_shape, T_inv_scale> poisson_gamma_lccdf(const T_n& n,
                                                      const T_shape& alpha,
                                                      const T_inv_scale& beta) {
  static const char* function = "poisson_gamma_lccdf";

  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_beta_ref = ref_type_t<T_inv_scale>;

  T_n_ref n_ref = n;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  check_nonnegative(function, "Random variable", n_ref);
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  return neg_binomial_2_lccdf<propto>(n_ref,
                                     elt_divide(alpha_ref, beta_ref),
                                     alpha_ref);
}

template <typename T_n, typename T_shape, typename T_inv_scale>
inline return_type_t<T_shape, T_inv_scale> poisson_gamma_lccdf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  return poisson_gamma_lccdf<false>(n, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
