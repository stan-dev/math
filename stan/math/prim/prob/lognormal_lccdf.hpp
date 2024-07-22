#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale> lognormal_lccdf(const T_y& y,
                                                   const T_loc& mu,
                                                   const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  static constexpr const char* function = "lognormal_lccdf";

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, mu, sigma)) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  if (sum(promote_scalar<int>(y_val == 0))) {
    return ops_partials.build(0.0);
  }

  const auto& log_y = log(y_val);
  const auto& scaled_diff
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
          (log_y - mu_val) / (sigma_val * SQRT_TWO));
  const auto& erfc_calc
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
          erfc(scaled_diff));

  size_t N = max_size(y, mu, sigma);
  T_partials_return ccdf_log = N * LOG_HALF + sum(log(erfc_calc));

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    const auto& exp_m_sq_diff = exp(-scaled_diff * scaled_diff);
    const auto& rep_deriv = to_ref_if<!is_constant_all<T_y>::value
                                          + !is_constant_all<T_scale>::value
                                          + !is_constant_all<T_loc>::value
                                      >= 2>(
        SQRT_TWO_OVER_SQRT_PI * exp_m_sq_diff / (sigma_val * erfc_calc));
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = -rep_deriv / y_val;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) = rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      partials<2>(ops_partials) = rep_deriv * scaled_diff * SQRT_TWO;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
