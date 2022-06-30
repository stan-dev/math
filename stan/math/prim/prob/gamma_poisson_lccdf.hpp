#ifndef STAN_MATH_PRIM_PROB_GAMMA_POISSON_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_GAMMA_POISSON_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/neg_binomial_2_lccdf.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_n, typename T_shape, typename T_inv_scale>
return_type_t<T_shape, T_inv_scale> gamma_poisson_lccdf(
    const T_n& n, const T_shape& alpha, const T_inv_scale& beta) {
  static const char* function = "gamma_poisson_lccdf";
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

  return neg_binomial_2_lccdf(n_ref, elt_divide(alpha_ref, beta_ref),
                              alpha_ref);
}

}  // namespace math
}  // namespace stan
#endif
