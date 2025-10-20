#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_cdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "skew_normal_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Shape paramter", alpha);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;
  T_alpha_ref alpha_ref = alpha;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", sigma_val);
  check_finite(function, "Shape parameter", alpha_val);

  if (size_zero(y, mu, sigma, alpha)) {
    return 1.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, alpha_ref);

  auto diff = to_ref((y_val - mu_val) / sigma_val);
  auto scaled_diff
      = to_ref_if<is_any_autodiff_v<T_y, T_loc, T_scale>>(diff / SQRT_TWO);
  auto erfc_m_scaled_diff = erfc(-scaled_diff);
  auto owens_t_diff_alpha = owens_t(diff, alpha_val);
  auto cdf_ = to_ref(0.5 * erfc_m_scaled_diff - 2 * owens_t_diff_alpha);

  T_partials_return cdf = prod(cdf_);

  if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale, T_shape>) {
    auto cdf_quot = cdf / cdf_;
    auto diff_square = square(diff);
    if constexpr (is_any_autodiff_v<T_y, T_loc, T_scale>) {
      auto erfc_m_alpha_scaled_diff = erfc(-alpha_val * scaled_diff);
      auto exp_m_scaled_diff_square = exp(-0.5 * diff_square);
      auto rep_deriv = erfc_m_alpha_scaled_diff * INV_SQRT_TWO_PI * cdf_quot
                       / sigma_val * std::move(exp_m_scaled_diff_square);
      if constexpr (is_autodiff_v<T_y>) {
        partials<0>(ops_partials) = rep_deriv;
      }
      if constexpr (is_autodiff_v<T_loc>) {
        partials<1>(ops_partials) = -rep_deriv;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        partials<2>(ops_partials) = -std::move(rep_deriv) * diff;
      }
    }
    if constexpr (is_autodiff_v<T_shape>) {
      auto alpha_square = square(alpha_val);
      auto exp_tmp = exp(-0.5 * std::move(diff_square) * (1.0 + alpha_square));
      edge<3>(ops_partials).partials_
          = -2.0 * std::move(exp_tmp) / ((1 + std::move(alpha_square)) * TWO_PI)
            * std::move(cdf_quot);
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
