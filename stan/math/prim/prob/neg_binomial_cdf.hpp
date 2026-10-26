#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/any.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inc_beta_ddz.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/prod.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <limits>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the CDF of the negative binomial distribution. Given containers of
 * matching sizes, returns the product of probabilities.
 *
 * With p = beta / (1 + beta), the CDF at n is the regularized incomplete beta
 * function I_p(alpha, n + 1).
 *
 * @tparam T_n type of integer parameter
 * @tparam T_shape type of shape parameter
 * @tparam T_inv_scale type of inverse scale parameter
 * @param n failures variable
 * @param alpha shape parameter
 * @param beta inverse scale parameter
 * @return probability or product of probabilities
 * @throw std::domain_error if alpha or beta is not positive and finite
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_n, typename T_shape, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_n, T_shape, T_inv_scale>* = nullptr>
inline return_type_t<T_shape, T_inv_scale> neg_binomial_cdf(
    T_n&& n, T_shape&& alpha, T_inv_scale&& beta) {
  using T_partials_return = partials_return_t<T_n, T_shape, T_inv_scale>;
  using T_n_ref = ref_type_t<T_n>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  using T_beta_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static constexpr const char* function = "neg_binomial_cdf";
  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  T_n_ref n_ref = std::forward<T_n>(n);
  T_alpha_ref alpha_ref = std::forward<T_shape>(alpha);
  T_beta_ref beta_ref = std::forward<T_inv_scale>(beta);

  decltype(auto) n_val = to_ref(as_value_column_array_or_scalar(n_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Inverse scale parameter", beta_val);

  if (size_zero(n_ref, alpha_ref, beta_ref)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(alpha_ref, beta_ref);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  if (any(n_val < 0)) {
    return ops_partials.build(0.0);
  }
  if (any(n_val == std::numeric_limits<int>::max())) {
    return ops_partials.build(1.0);
  }

  if constexpr (is_autodiff_v<T_shape> || is_autodiff_v<T_inv_scale>) {
    auto inv_beta_p1 = to_ref(inv(beta_val + 1.0));
    auto p_val = to_ref(beta_val * inv_beta_p1);
    auto n_p1 = to_ref(promote_scalar<T_partials_return>(n_val + 1.0));
    auto P_i = to_ref(as_array_or_scalar(inc_beta(alpha_val, n_p1, p_val)));
    const T_partials_return P = prod(P_i);
    // d/dx prod(P_i) = sum(dP_i/dx * P / P_i)
    if constexpr (is_autodiff_v<T_shape> && is_autodiff_v<T_inv_scale>) {
      auto P_div_P_i = to_ref(P / P_i);
      partials<0>(ops_partials)
          = inc_beta_dda(alpha_val, n_p1, p_val, to_ref(digamma(alpha_val)),
                         to_ref(digamma(alpha_val + n_p1)))
            * P_div_P_i;
      partials<1>(ops_partials)
          = as_array_or_scalar(inc_beta_ddz(alpha_val, n_p1, p_val))
            * square(inv_beta_p1) * P_div_P_i;
    } else if constexpr (is_autodiff_v<T_shape>) {
      partials<0>(ops_partials)
          = inc_beta_dda(alpha_val, n_p1, p_val, to_ref(digamma(alpha_val)),
                         to_ref(digamma(alpha_val + n_p1)))
            / P_i * P;
    } else if constexpr (is_autodiff_v<T_inv_scale>) {
      partials<1>(ops_partials)
          = as_array_or_scalar(inc_beta_ddz(alpha_val, n_p1, p_val))
            * square(inv_beta_p1) / P_i * P;
    }
    return ops_partials.build(P);
  } else {
    auto inv_beta_p1 = to_ref(inv(beta_val + 1.0));
    auto p_val = to_ref(beta_val * inv_beta_p1);
    auto n_p1 = to_ref(promote_scalar<T_partials_return>(n_val + 1.0));
    const T_partials_return P
        = prod(as_array_or_scalar(inc_beta(alpha_val, n_p1, p_val)));
    return ops_partials.build(P);
  }
}

}  // namespace math
}  // namespace stan
#endif
