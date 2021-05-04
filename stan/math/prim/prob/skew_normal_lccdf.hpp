#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/owens_t.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  static const char* function = "skew_normal_lccdf";
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
    return 0;
  }

  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref, T_alpha_ref>
      ops_partials(y_ref, mu_ref, sigma_ref, alpha_ref);

  const auto& diff = to_ref((y_val - mu_val) / sigma_val);
  const auto& scaled_diff
      = to_ref_if<!is_constant_all<T_y, T_scale, T_loc>::value>(diff
                                                                / SQRT_TWO);
  const auto& erfc_m_scaled_diff = erfc(-scaled_diff);
  const auto& owens_t_diff_alpha = owens_t(diff, alpha_val);
  const auto& ccdf_log_ = to_ref_if<!is_constant_all<T_shape>::value>(
      1.0 - 0.5 * erfc_m_scaled_diff + 2 * owens_t_diff_alpha);

  T_partials_return ccdf_log = sum(log(ccdf_log_));

  if (!is_constant_all<T_y, T_scale, T_loc, T_shape>::value) {
    const auto& diff_square
        = to_ref_if<(!is_constant_all<T_y, T_scale, T_loc>::value
                     && !is_constant_all<T_shape>::value)>(square(diff));
    if (!is_constant_all<T_y, T_scale, T_loc>::value) {
      const auto& erf_alpha_scaled_diff = erf(alpha_val * scaled_diff);
      const auto& exp_m_scaled_diff_square = exp(-0.5 * diff_square);
      auto rep_deriv = to_ref_if<!is_constant_all<T_y>::value
                                     + !is_constant_all<T_scale>::value
                                     + !is_constant_all<T_loc>::value
                                 >= 2>(
          (erf_alpha_scaled_diff + 1) * INV_SQRT_TWO_PI
          / (ccdf_log_ * sigma_val) * exp_m_scaled_diff_square);
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_ = -rep_deriv;
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge3_.partials_ = rep_deriv * diff;
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge2_.partials_ = std::move(rep_deriv);
      }
    }
    if (!is_constant_all<T_shape>::value) {
      const auto& alpha_square = square(alpha_val);
      ops_partials.edge4_.partials_
          = 2.0 * exp(-0.5 * diff_square * (1.0 + alpha_square))
            / ((1 + alpha_square) * TWO_PI * ccdf_log_);
    }
  }

  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
