#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
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

// Loglogistic(y | alpha, beta)    [y > 0, alpha > 0, beta > 0]
template <bool propto, typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> loglogistic_lpdf(const T_y& y,
                                                      const T_scale& alpha,
                                                      const T_shape& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_scale_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using T_shape_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using std::pow;
  static const char* function = "loglogistic_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         alpha, "Shape parameter", beta);

  T_y_ref y_ref = y;
  T_scale_ref alpha_ref = alpha;
  T_shape_ref beta_ref = beta;

  // Tukaj moram uporabit ta partials_return_t? Namesto auto? Hm...
  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_positive_finite(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", alpha_val);
  check_positive_finite(function, "Shape parameter", beta_val);

  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }
  if (!include_summand<propto, T_y, T_scale, T_shape>::value) {
    return 0.0;
  }

  operands_and_partials<T_y_ref, T_scale_ref, T_shape_ref> ops_partials(
      y_ref, alpha_ref, beta_ref);

  const auto& inv_alpha
      = to_ref_if<!is_constant_all<T_y, T_scale>::value>(inv(alpha_val));
  const auto& y_div_alpha = to_ref(y_val * inv_alpha);
  const auto& y_div_alpha_pow_beta =
      to_ref_if<!is_constant_all<T_shape>::value>(pow(y_div_alpha, beta_val));
  const auto& log1_arg =
      to_ref_if<!is_constant_all<T_y, T_scale, T_shape>::value>(1 +
      y_div_alpha_pow_beta);
  const auto& log_y = to_ref_if<!is_constant_all<T_shape>::value>(log(y_val));

  size_t N = max_size(y, alpha, beta);
  size_t N_alpha_beta = max_size(alpha, beta);

  T_partials_return logp = sum((beta_val - 1) * log_y - 2 * log1_arg);
  // T_partials_return logp = sum((beta_val - 1) * log(y_val) -
  //   2 * log1p(pow((y_val * inv(alpha_val)), beta_val)));

  if (include_summand<propto, T_scale, T_shape>::value) {
    logp += sum(N * (log(beta_val) - log(alpha_val) -
      (beta_val - 1) * log(alpha_val)) / N_alpha_beta);
  }

  if (!is_constant_all<T_y, T_scale, T_shape>::value) {
    const auto& two_inv_log1_arg = to_ref_if<!is_constant_all<T_y>::value
                                             + !is_constant_all<T_scale>::value
                                             + !is_constant_all<T_shape>::value
                                             >= 2>(2 * inv(log1_arg));
    const auto& y_pow_beta =
        to_ref_if<!is_constant_all<T_y, T_scale>::value>(pow(y_val, beta_val));

    if (!is_constant_all<T_y>::value) {
      const auto& y_deriv = (beta_val - 1.0) * inv(y_val) -
        two_inv_log1_arg *
        (beta_val * inv(pow(alpha_val, beta_val))) * y_pow_beta * inv(y_val);
        ops_partials.edge1_.partials_ = y_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      const auto& alpha_deriv = - beta_val * inv(alpha_val) -
        two_inv_log1_arg *
        y_pow_beta * (-beta_val) * pow(alpha_val, -beta_val - 1);
      ops_partials.edge2_.partials_ = alpha_deriv;
    }
    if (!is_constant_all<T_shape>::value) {
      const auto& beta_deriv = (1.0 * inv(beta_val)) + log(y_val) - log(alpha_val) -
        two_inv_log1_arg *
        pow((y_val * inv(alpha_val)), beta_val) * log(y_val * inv(alpha_val));
      ops_partials.edge3_.partials_ = beta_deriv;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_scale, T_shape> loglogistic_lpdf(const T_y& y,
                                                        const T_scale& alpha,
                                                        const T_shape& beta) {
  return loglogistic_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
