#ifndef STAN_MATH_PRIM_PROB_PARETO_LCDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_LCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_scale, typename T_shape>
return_type_t<T_y, T_scale, T_shape> pareto_lcdf(const T_y& y,
                                                 const T_scale& y_min,
                                                 const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_y_min_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using std::isinf;
  static const char* function = "pareto_lcdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         y_min, "Shape parameter", alpha);

  if (size_zero(y, y_min, alpha)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_y_min_ref y_min_ref = y_min;
  T_alpha_ref alpha_ref = alpha;

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& y_min_col = as_column_vector_or_scalar(y_min_ref);
  const auto& alpha_col = as_column_vector_or_scalar(alpha_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& y_min_arr = as_array_or_scalar(y_min_col);
  const auto& alpha_arr = as_array_or_scalar(alpha_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(y_min_arr))> y_min_val = value_of(y_min_arr);
  ref_type_t<decltype(value_of(alpha_arr))> alpha_val = value_of(alpha_arr);

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", y_min_val);
  check_positive_finite(function, "Shape parameter", alpha_val);

  operands_and_partials<T_y_ref, T_y_min_ref, T_alpha_ref> ops_partials(
      y_ref, y_min_ref, alpha_ref);

  // Explicit return for extreme values
  // The gradients are technically ill-defined, but treated as zero
  if (sum(promote_scalar<int>(y_val < y_min_val))) {
    return ops_partials.build(negative_infinity());
  }
  if (sum(promote_scalar<int>(isinf(y_val)))) {
    return ops_partials.build(0.0);
  }

  const auto& log_quot
      = to_ref_if<!is_constant_all<T_y, T_scale, T_shape>::value>(
          log(y_min_val / y_val));
  const auto& exp_prod
      = to_ref_if<!is_constant_all<T_y, T_scale, T_shape>::value>(
          exp(alpha_val * log_quot));
  T_partials_return P = sum(log(1 - exp_prod));

  if (!is_constant_all<T_y, T_scale, T_shape>::value) {
    const auto& common_deriv = to_ref_if<(!is_constant_all<T_y, T_scale>::value
                                          && !is_constant_all<T_shape>::value)>(
        exp_prod / (1 - exp_prod));
    if (!is_constant_all<T_y, T_scale>::value) {
      const auto& y_min_inv = inv(y_min_val);
      const auto& common_deriv2 = to_ref_if<(
          !is_constant_all<T_y>::value && !is_constant_all<T_scale>::value)>(
          -alpha_val * y_min_inv * common_deriv);
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_ = -common_deriv2 * exp(log_quot);
      }
      if (!is_constant_all<T_scale>::value) {
        ops_partials.edge2_.partials_ = std::move(common_deriv2);
      }
    }
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge3_.partials_ = -common_deriv * log_quot;
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
