#ifndef STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_PARETO_TYPE_2_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/pow.hpp>
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
return_type_t<T_y, T_loc, T_scale, T_shape> pareto_type_2_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& lambda,
    const T_shape& alpha) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_mu_ref = ref_type_if_not_constant_t<T_loc>;
  using T_lambda_ref = ref_type_if_not_constant_t<T_scale>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_shape>;
  static constexpr const char* function = "pareto_type_2_lccdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", lambda, "Shape parameter",
                         alpha);

  if (size_zero(y, mu, lambda, alpha)) {
    return 0;
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

  auto ops_partials
      = make_partials_propagator(y_ref, mu_ref, lambda_ref, alpha_ref);

  const auto& log_temp = to_ref_if<!is_constant_all<T_shape>::value>(
      log1p((y_val - mu_val) / lambda_val));
  T_partials_return P = -sum(alpha_val * log_temp);

  if (!is_constant_all<T_y, T_loc, T_scale>::value) {
    auto rep_deriv = to_ref_if<(!is_constant_all<T_y>::value
                                + !is_constant_all<T_scale>::value
                                + !is_constant_all<T_loc>::value)
                               >= 2>(alpha_val / (y_val - mu_val + lambda_val));
    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) = -rep_deriv;
    }
    if (!is_constant_all<T_scale>::value) {
      edge<2>(ops_partials).partials_
          = rep_deriv * (y_val - mu_val) / lambda_val;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) = std::move(rep_deriv);
    }
  }
  size_t N = max_size(y, mu, lambda, alpha);
  if (!is_constant_all<T_shape>::value) {
    if (is_vector<T_shape>::value) {
      using Log_temp_scalar = partials_return_t<T_y, T_loc, T_scale>;
      using Log_temp_array = Eigen::Array<Log_temp_scalar, Eigen::Dynamic, 1>;
      if (is_vector<T_y>::value || is_vector<T_loc>::value
          || is_vector<T_scale>::value) {
        partials<3>(ops_partials) = -forward_as<Log_temp_array>(log_temp);
      } else {
        partials<3>(ops_partials) = Log_temp_array::Constant(
            N, 1, -forward_as<Log_temp_scalar>(log_temp));
      }
    } else {
      forward_as<internal::broadcast_array<T_partials_return>>(
          partials<3>(ops_partials))
          = -log_temp * N / max_size(y, mu, lambda);
    }
  }

  return ops_partials.build(P);
}

}  // namespace math
}  // namespace stan
#endif
