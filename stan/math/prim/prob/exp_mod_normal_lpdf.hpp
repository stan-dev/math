#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static const char* function = "exp_mod_normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale paramter",
                         lambda);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_lambda_ref lambda_ref = lambda;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);
  check_positive_finite(function, "Inv_scale parameter", lambda_val);

  if (size_zero(y, mu, sigma, lambda)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_loc, T_scale, T_inv_scale>::value) {
    return 0.0;
  }

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_loc>::value>(inv(sigma_val));
  const auto& sigma_sq
      = to_ref_if<!is_constant_all<T_scale>::value>(square(sigma_val));
  const auto& lambda_sigma_sq = to_ref(lambda_val * sigma_sq);
  const auto& mu_minus_y = to_ref(mu_val - y_val);
  const auto& inner_term
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale, T_inv_scale>::value>(
          (mu_minus_y + lambda_sigma_sq) * INV_SQRT_TWO * inv_sigma);
  const auto& erfc_calc = to_ref(erfc(inner_term));

  size_t N = max_size(y, mu, sigma, lambda);
  T_partials_return logp(0.0);
  if (include_summand<propto>::value) {
    logp -= LOG_TWO * N;
  }
  if (include_summand<propto, T_inv_scale>::value) {
    logp += sum(log(lambda_val)) * N / math::size(lambda);
  }
  const auto& log_erfc_calc = log(erfc_calc);
  logp
      += sum(lambda_val * (mu_minus_y + 0.5 * lambda_sigma_sq) + log_erfc_calc);

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);

  if (!is_constant_all<T_y, T_loc, T_scale, T_inv_scale>::value) {
    const auto& exp_m_sq_inner_term = exp(-square(inner_term));
    const auto& deriv_logerfc = to_ref_if<
        !is_constant_all<T_y, T_loc>::value + !is_constant_all<T_scale>::value
            + !is_constant_all<T_inv_scale>::value
        >= 2>(-SQRT_TWO_OVER_SQRT_PI * exp_m_sq_inner_term / erfc_calc);
    if (!is_constant_all<T_y, T_loc>::value) {
      const auto& deriv = to_ref_if < !is_constant_all<T_y>::value
                          && !is_constant_all<T_loc>::value
                                 > (lambda_val + deriv_logerfc * inv_sigma);
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials) = -deriv;
      }
      if (!is_constant_all<T_loc>::value) {
        partials<1>(ops_partials) = deriv;
      }
    }
    if (!is_constant_all<T_scale>::value) {
      edge<2>(ops_partials).partials_
          = sigma_val * square(lambda_val)
            + deriv_logerfc * (lambda_val - mu_minus_y / sigma_sq);
    }
    if (!is_constant_all<T_inv_scale>::value) {
      partials<3>(ops_partials) = inv(lambda_val) + lambda_sigma_sq + mu_minus_y
                                  + deriv_logerfc * sigma_val;
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  return exp_mod_normal_lpdf<false>(y, mu, sigma, lambda);
}

}  // namespace math
}  // namespace stan
#endif
