#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LPDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

// pareto_type_2(y|lambda, alpha)  [y >= 0;  lambda > 0;  alpha > 0]
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_lambda_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  static const char* function = "pareto_type_2_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", lambda, "Shape parameter",
                         alpha);

  if (size_zero(y, mu, lambda, alpha)) {
    return 0.0;
  }

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_lambda_ref lambda_ref = lambda;
  T_alpha_ref alpha_ref = alpha;

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& mu_col = as_column_vector_or_scalar(mu_ref);
  const auto& lambda_col = as_column_vector_or_scalar(lambda_ref);
  const auto& alpha_col = as_column_vector_or_scalar(alpha_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& mu_arr = as_array_or_scalar(mu_col);
  const auto& lambda_arr = as_array_or_scalar(lambda_col);
  const auto& alpha_arr = as_array_or_scalar(alpha_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(mu_arr))> mu_val = value_of(mu_arr);
  ref_type_t<decltype(value_of(lambda_arr))> lambda_val = value_of(lambda_arr);
  ref_type_t<decltype(value_of(alpha_arr))> alpha_val = value_of(alpha_arr);

  check_greater_or_equal(function, "Random variable", y_val, mu_val);
  check_positive_finite(function, "Scale parameter", lambda_val);
  check_positive_finite(function, "Shape parameter", alpha_val);

  if (!include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    return 0.0;
  }

  const auto& log1p_scaled_diff = to_ref_if<!is_constant_all<T_shape>::value>(
      log1p((y_val - mu_val) / lambda_val));

  size_t N = max_size(y, mu, lambda, alpha);
  T_partials_return logp(0.0);
  if (include_summand<propto, T_shape>::value) {
    logp += sum(log(alpha_val)) * N / math::size(alpha);
  }
  if (include_summand<propto, T_scale>::value) {
    logp -= sum(log(lambda_val)) * N / math::size(lambda);
  }
  if (include_summand<propto, T_y, T_loc, T_scale, T_shape>::value) {
    logp -= sum((alpha_val + 1.0) * log1p_scaled_diff);
  }

  operands_and_partials<T_y_ref, T_mu_ref, T_lambda_ref, T_alpha_ref>
      ops_partials(y_ref, mu_ref, lambda_ref, alpha_ref);

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    const auto& inv_sum = to_ref_if<(!is_constant_all<T_y, T_loc>::value
                                     && !is_constant_all<T_scale>::value)>(
        inv(lambda_val + y_val - mu_val));
    const auto& alpha_div_sum
        = to_ref_if<(!is_constant_all<T_y, T_loc>::value
                     && !is_constant_all<T_scale>::value)>(alpha_val * inv_sum);
    if (!is_constant_all<T_y, T_loc>::value) {
      const auto& deriv_1_2 = to_ref_if<(!is_constant_all<T_y>::value
                                         && !is_constant_all<T_loc>::value)>(
          inv_sum + alpha_div_sum);
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_ = -deriv_1_2;
      }
      if (!is_constant_all<T_loc>::value) {
        ops_partials.edge2_.partials_ = std::move(deriv_1_2);
      }
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_
          = alpha_div_sum * (y_val - mu_val) / lambda_val - inv_sum;
    }
  }
  if (!is_constant_all<T_shape>::value) {
    ops_partials.edge4_.partials_ = inv(alpha_val) - log1p_scaled_diff;
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  return pareto_type_2_lpdf<false>(y, mu, lambda, alpha);
}

}  // namespace math
}  // namespace stan
#endif
