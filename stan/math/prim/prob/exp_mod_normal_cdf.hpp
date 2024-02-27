#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_inv_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_cdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;
  using T_partials_array = Eigen::Array<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_sigma_ref = ref_type_if_not_constant_t<T_scale>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_inv_scale>;
  static constexpr const char* function = "exp_mod_normal_cdf";
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
    return 1.0;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, sigma_ref, lambda_ref);

  using T_y_val_scalar = scalar_type_t<decltype(y_val)>;
  if (is_vector<T_y>::value) {
    if ((forward_as<Eigen::Array<T_y_val_scalar, Eigen::Dynamic, 1>>(y_val)
         == NEGATIVE_INFTY)
            .any()) {
      return ops_partials.build(0.0);
    }
  } else {
    if (forward_as<T_y_val_scalar>(y_val) == NEGATIVE_INFTY) {
      return ops_partials.build(0.0);
    }
  }

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(inv(sigma_val));
  const auto& diff = to_ref(y_val - mu_val);
  const auto& v = to_ref(lambda_val * sigma_val);
  const auto& scaled_diff = to_ref(diff * INV_SQRT_TWO * inv_sigma);
  const auto& scaled_diff_diff
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale, T_inv_scale>::value>(
          scaled_diff - v * INV_SQRT_TWO);
  const auto& erf_calc = to_ref(0.5 * (1 + erf(scaled_diff_diff)));

  const auto& exp_term
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale, T_inv_scale>::value>(
          exp(0.5 * square(v) - lambda_val * diff));
  const auto& cdf_n
      = to_ref(0.5 + 0.5 * erf(scaled_diff) - exp_term * erf_calc);

  T_partials_return cdf(1.0);
  if (is_vector<decltype(cdf_n)>::value) {
    cdf = forward_as<T_partials_array>(cdf_n).prod();
  } else {
    cdf = forward_as<T_partials_return>(cdf_n);
  }

  if (!is_constant_all<T_y, T_loc, T_scale, T_inv_scale>::value) {
    const auto& exp_term_2
        = to_ref_if<(!is_constant_all<T_y, T_loc, T_scale>::value
                     && !is_constant_all<T_inv_scale>::value)>(
            exp(-square(scaled_diff_diff)));
    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      constexpr bool need_deriv_refs = !is_constant_all<T_y, T_loc>::value
                                       && !is_constant_all<T_scale>::value;
      const auto& deriv_1
          = to_ref_if<need_deriv_refs>(lambda_val * exp_term * erf_calc);
      const auto& deriv_2 = to_ref_if<need_deriv_refs>(
          INV_SQRT_TWO_PI * exp_term * exp_term_2 * inv_sigma);
      const auto& sq_scaled_diff = square(scaled_diff);
      const auto& exp_m_sq_scaled_diff = exp(-sq_scaled_diff);
      const auto& deriv_3 = to_ref_if<need_deriv_refs>(
          INV_SQRT_TWO_PI * exp_m_sq_scaled_diff * inv_sigma);
      if (!is_constant_all<T_y, T_loc>::value) {
        const auto& deriv = to_ref_if<(!is_constant_all<T_loc>::value
                                       && !is_constant_all<T_y>::value)>(
            cdf * (deriv_1 - deriv_2 + deriv_3) / cdf_n);
        if (!is_constant_all<T_y>::value) {
          partials<0>(ops_partials) = deriv;
        }
        if (!is_constant_all<T_loc>::value) {
          partials<1>(ops_partials) = -deriv;
        }
      }
      if (!is_constant_all<T_scale>::value) {
        edge<2>(ops_partials).partials_
            = -cdf
              * ((deriv_1 - deriv_2) * v
                 + (deriv_3 - deriv_2) * scaled_diff * SQRT_TWO)
              / cdf_n;
      }
    }
    if (!is_constant_all<T_inv_scale>::value) {
      edge<3>(ops_partials).partials_
          = cdf * exp_term
            * (INV_SQRT_TWO_PI * sigma_val * exp_term_2
               - (v * sigma_val - diff) * erf_calc)
            / cdf_n;
    }
  }
  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
