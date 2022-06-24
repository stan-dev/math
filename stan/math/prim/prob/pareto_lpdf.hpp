#ifndef STAN_MATH_PRIM_PROB_PARETO_LPDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
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

// Pareto(y|y_m, alpha)  [y > y_m;  y_m > 0;  alpha > 0]
template <bool propto, typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> pareto_lpdf(const T_y& y,
                                                 const T_scale& y_min,
                                                 const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_y_min_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  static const char* function = "pareto_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         y_min, "Shape parameter", alpha);
  if (size_zero(y, y_min, alpha)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_y_min_ref y_min_ref = y_min;
  T_alpha_ref alpha_ref = alpha;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) y_min_val = to_ref(as_value_column_array_or_scalar(y_min_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));

  check_not_nan(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", y_min_val);
  check_positive_finite(function, "Shape parameter", alpha_val);

  if (!include_summand<propto, T_y, T_scale, T_shape>::value) {
    return 0;
  }

  if (sum(promote_scalar<int>(y_val < y_min_val))) {
    return LOG_ZERO;
  }

  const auto& log_y = to_ref_if<!is_constant_all<T_shape>::value>(log(y_val));

  size_t N = max_size(y, y_min, alpha);
  T_partials_return logp(0);
  if (include_summand<propto, T_shape>::value) {
    logp = sum(log(alpha_val)) * N / math::size(alpha);
  }
  if (include_summand<propto, T_y, T_shape>::value) {
    logp -= sum(alpha_val * log_y + log_y) * N / max_size(alpha, y);
  }

  operands_and_partials<T_y_ref, T_y_min_ref, T_alpha_ref> ops_partials(
      y_ref, y_min_ref, alpha_ref);
  if (!is_constant_all<T_y>::value) {
    const auto& inv_y = inv(y_val);
    ops_partials.edge1_.partials_
        = -(alpha_val * inv_y + inv_y) * N / max_size(alpha, y);
  }
  if (!is_constant_all<T_scale>::value) {
    ops_partials.edge2_.partials_
        = alpha_val / y_min_val * N / max_size(alpha, y_min);
  }
  if (include_summand<propto, T_scale, T_shape>::value) {
    const auto& log_y_min
        = to_ref_if<!is_constant_all<T_shape>::value>(log(y_min_val));
    logp += sum(alpha_val * log_y_min) * N / max_size(alpha, y_min);
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge3_.partials_ = inv(alpha_val) + log_y_min - log_y;
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_scale, T_shape> pareto_lpdf(const T_y& y,
                                                        const T_scale& y_min,
                                                        const T_shape& alpha) {
  return pareto_lpdf<false>(y, y_min, alpha);
}

}  // namespace math
}  // namespace stan
#endif
