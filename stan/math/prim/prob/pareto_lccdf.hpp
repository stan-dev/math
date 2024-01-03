#ifndef STAN_MATH_PRIM_PROB_PARETO_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
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

template <typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> pareto_lccdf(const T_y& y,
                                                  const T_scale& y_min,
                                                  const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_y_min_ref = ref_type_if_not_constant_t<T_scale>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  using std::isinf;
  static const char* function = "pareto_lccdf";
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

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", y_min_val);
  check_positive_finite(function, "Shape parameter", alpha_val);

  auto ops_partials = make_partials_propagator(y_ref, y_min_ref, alpha_ref);

  if (sum(promote_scalar<int>(y_val < y_min_val))) {
    return ops_partials.build(0.0);
  }
  if (sum(promote_scalar<int>(isinf(y_val)))) {
    return ops_partials.build(negative_infinity());
  }

  auto log_quot = to_ref_if<(!is_constant_all<T_y>::value
                             || !is_constant_all<T_shape>::value)>(
      log(y_min_val / y_val));

  T_partials_return P = sum(alpha_val * log_quot);

  size_t N = max_size(y, y_min, alpha);
  if (!is_constant_all<T_y, T_scale>::value) {
    const auto& alpha_div_y_min = to_ref_if<(
        !is_constant_all<T_y>::value && !is_constant_all<T_scale>::value)>(
        alpha_val / y_min_val);
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = -alpha_div_y_min * exp(log_quot);
    }
    if (!is_constant_all<T_scale>::value) {
      edge<1>(ops_partials).partials_
          = alpha_div_y_min * N / max_size(y_min, alpha);
    }
  }
  if (!is_constant_all<T_shape>::value) {
    if (is_vector<T_shape>::value) {
      using Log_quot_scalar = partials_return_t<T_y, T_scale>;
      using Log_quot_array = Eigen::Array<Log_quot_scalar, Eigen::Dynamic, 1>;
      if (is_vector<T_y>::value || is_vector<T_scale>::value) {
        edge<2>(ops_partials).partials_
            = forward_as<Log_quot_array>(std::move(log_quot));
      } else {
        partials<2>(ops_partials) = Log_quot_array::Constant(
            N, 1, forward_as<Log_quot_scalar>(log_quot));
      }
    } else {
      forward_as<internal::broadcast_array<T_partials_return>>(
          partials<2>(ops_partials))
          = log_quot * N / max_size(y, y_min);
    }
  }
  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
