#ifndef STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGNORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/prod.hpp>
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
return_type_t<T_y, T_loc, T_scale> lognormal_cdf(const T_y& y, const T_loc& mu,
                                                 const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "lognormal_cdf";

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
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, sigma_ref);

  if (sum(promote_scalar<int>(y_val == 0))) {
    return ops_partials.build(0.0);
  }

  const auto& log_y = log(y_val);
  const auto& scaled_diff
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
          (log_y - mu_val) / (sigma_val * SQRT_TWO));
  const auto& erfc_m_diff = erfc(-scaled_diff);
  const auto& cdf_n = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
      0.5 * erfc_m_diff);

  T_partials_return cdf = prod(cdf_n);

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    const auto& exp_m_sq_diff = exp(-scaled_diff * scaled_diff);
    const auto& rep_deriv = to_ref_if<!is_constant_all<T_y>::value
                                          + !is_constant_all<T_scale>::value
                                          + !is_constant_all<T_loc>::value
                                      >= 2>(
        -cdf * INV_SQRT_TWO_PI * exp_m_sq_diff / (sigma_val * cdf_n));
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
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
