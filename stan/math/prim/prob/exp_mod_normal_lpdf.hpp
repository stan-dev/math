#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/std_normal_lcdf_impl.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lpdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_inv_scale&& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static constexpr const char* function = "exp_mod_normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale parameter",
                         lambda);
  T_y_ref y_ref = std::forward<T_y>(y);
  T_mu_ref mu_ref = std::forward<T_loc>(mu);
  T_sigma_ref sigma_ref = std::forward<T_scale>(sigma);
  T_lambda_ref lambda_ref = std::forward<T_inv_scale>(lambda);

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);
  check_positive_finite(function, "Inv_scale parameter", lambda_val);

  if (size_zero(y_ref, mu_ref, sigma_ref, lambda_ref)) {
    return 0.0;
  }
  if constexpr (!include_summand<propto, T_y, T_loc, T_scale,
                                 T_inv_scale>::value) {
    return 0.0;
  }

  const auto& inv_sigma
      = to_ref_if<is_any_autodiff_v<T_y, T_loc>>(inv(sigma_val));
  const auto& sigma_sq = to_ref_if<is_autodiff_v<T_scale>>(square(sigma_val));
  const auto& lambda_sigma_sq = to_ref(lambda_val * sigma_sq);
  const auto& mu_minus_y = to_ref(mu_val - y_val);
  // log(erfc(t) / 2) = log Phi(-sqrt(2) t) cancels the log(1/2) constant.
  const auto& z = to_ref(-(mu_minus_y + lambda_sigma_sq) * inv_sigma);
  const auto [values, slopes] = internal::std_normal_lcdf_value_grad<
      is_any_autodiff_v<T_y, T_loc, T_scale, T_inv_scale>>(z);

  size_t N = max_size(y_ref, mu_ref, sigma_ref, lambda_ref);
  T_partials_return logp(0.0);
  if constexpr (include_summand<propto, T_inv_scale>::value) {
    logp += sum(log(lambda_val)) * N / math::size(lambda_ref);
  }
  logp += sum(lambda_val * (mu_minus_y + 0.5 * lambda_sigma_sq) + values);

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);

  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale, T_inv_scale>) {
    if constexpr (is_any_autodiff_v<T_y, T_loc>) {
      const auto& deriv = to_ref_if<is_all_autodiff_v<T_y, T_loc>>(
          lambda_val - slopes * inv_sigma);
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials) = -deriv;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials) = deriv;
      }
    }
    if constexpr (is_autodiff_v<T_scale>) {
      edge<2>(ops_partials).partials_
          = sigma_val * square(lambda_val)
            - slopes * (lambda_val - mu_minus_y * square(inv_sigma));
    }
    if constexpr (is_autodiff_v<T_inv_scale>) {
      partials<3>(ops_partials)
          = inv(lambda_val) + lambda_sigma_sq + mu_minus_y - slopes * sigma_val;
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lpdf(
    T_y&& y, T_loc&& mu, T_scale&& sigma, T_inv_scale&& lambda) {
  return exp_mod_normal_lpdf<false>(
      std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma), std::forward<T_inv_scale>(lambda));
}

}  // namespace math
}  // namespace stan
#endif
