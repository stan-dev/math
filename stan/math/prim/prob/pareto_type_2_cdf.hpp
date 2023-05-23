#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_CDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_loc, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_cdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_lambda_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using std::pow;
  static const char* function = "pareto_type_2_cdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", lambda, "Shape parameter",
                         alpha);

  if (size_zero(y, mu, lambda, alpha)) {
    return 1.0;
  }

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_lambda_ref lambda_ref = lambda;
  T_alpha_ref alpha_ref = alpha;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) lambda_val
      = to_ref(as_value_column_array_or_scalar(lambda_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", lambda_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_greater_or_equal(function, "Random variable", y_val, mu_val);

  const auto& summed = to_ref_if<!is_constant_all<T_y, T_loc, T_scale>::value>(
      lambda_val + y_val - mu_val);
  const auto& temp
      = to_ref_if<!is_constant_all<T_shape>::value>(summed / lambda_val);
  const auto& p1_pow_alpha
      = to_ref_if<!is_constant_all<T_y, T_loc, T_scale, T_shape>::value>(
          pow(temp, -alpha_val));
  T_partials_return P = prod(1.0 - p1_pow_alpha);

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, lambda_ref, alpha_ref);

  if (!is_constant_all<T_y, T_loc, T_scale, T_shape>::value) {
    const auto& P_div_Pn
        = to_ref_if<(!is_constant_all<T_y, T_loc, T_scale>::value
                     && !is_constant_all<T_shape>::value)>(
            P / (1.0 - p1_pow_alpha));
    if (!is_constant_all<T_y, T_loc, T_scale>::value) {
      auto grad_1_2
          = to_ref_if<!is_constant_all<T_loc>::value
                          + !is_constant_all<T_scale>::value
                          + !is_constant_all<T_y>::value
                      >= 2>(p1_pow_alpha / summed * alpha_val * P_div_Pn);
      if (!is_constant_all<T_loc>::value) {
        partials<1>(ops_partials) = -grad_1_2;
      }
      if (!is_constant_all<T_scale>::value) {
        edge<2>(ops_partials).partials_
            = (mu_val - y_val) / lambda_val * grad_1_2;
      }
      if (!is_constant_all<T_y>::value) {
        partials<0>(ops_partials) = std::move(grad_1_2);
      }
    }
    if (!is_constant_all<T_shape>::value) {
      partials<3>(ops_partials) = log(temp) * p1_pow_alpha * P_div_Pn;
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
