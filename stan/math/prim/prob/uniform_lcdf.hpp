#ifndef STAN_MATH_PRIM_PROB_UNIFORM_LCDF_HPP
#define STAN_MATH_PRIM_PROB_UNIFORM_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_low, typename T_high,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_low, T_high>* = nullptr>
return_type_t<T_y, T_low, T_high> uniform_lcdf(const T_y& y, const T_low& alpha,
                                               const T_high& beta) {
  using T_partials_return = partials_return_t<T_y, T_low, T_high>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_low>::value, T_low>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_high>::value, T_high>;
  static const char* function = "uniform_lcdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Lower bound parameter", alpha,
                         "Upper bound parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Lower bound parameter", alpha_val);
  check_finite(function, "Upper bound parameter", beta_val);
  check_greater(function, "Upper bound parameter", beta_val, alpha_val);

  if (size_zero(y, alpha, beta)) {
    return 0;
  }

  if (sum(promote_scalar<int>(y_val < alpha_val))
      || sum(promote_scalar<int>(beta_val < y_val))) {
    return negative_infinity();
  }

  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(
      y_ref, alpha_ref, beta_ref);

  const auto& b_minus_a
      = to_ref_if<!is_constant_all<T_y, T_low, T_high>::value>(beta_val
                                                               - alpha_val);
  const auto& y_minus_alpha
      = to_ref_if<!is_constant_all<T_y, T_low>::value>(y_val - alpha_val);
  const auto& cdf_log_n = y_minus_alpha / b_minus_a;
  T_partials_return cdf_log = sum(log(cdf_log_n));

  if (!is_constant_all<T_y>::value) {
    if (!is_vector<T_y>::value && is_vector<T_high>::value
        && !is_vector<T_low>::value) {
      ops_partials.edge1_.partials_ = size(beta) * inv(y_minus_alpha);
    } else {
      ops_partials.edge1_.partials_ = inv(y_minus_alpha);
    }
  }
  if (!is_constant_all<T_low>::value) {
    ops_partials.edge2_.partials_
        = (y_val - beta_val) / (b_minus_a * y_minus_alpha);
  }
  if (!is_constant_all<T_high>::value) {
    if (is_vector<T_y>::value && !is_vector<T_low>::value
        && !is_vector<T_high>::value) {
      ops_partials.edge3_.partials_ = inv(-b_minus_a) * size(y);
    } else {
      ops_partials.edge3_.partials_ = inv(-b_minus_a);
    }
  }
  return ops_partials.build(cdf_log);
}

}  // namespace math
}  // namespace stan
#endif
