#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_CDF_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <cmath>
#include <iostream>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The loglogistic cumulative distribution function for the specified
 * scalar(s) given the specified scales(s) and shape(s). y, alpha, or
 * beta can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 *
 * @tparam T_y type of scalar.
 * @tparam T_scale type of scale parameter.
 * @tparam T_shape type of shape parameter.
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) scale parameter(s)
 * for the loglogistic distribution.
 * @param beta (Sequence of) shape parameter(s) for the
 * loglogistic distribution.
 * @return The loglogistic cdf evaluated at the specified arguments.
 * @throw std::domain_error if any of the inputs are not positive or
 * if and of the parameters are not finite.
 */
template <typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
return_type_t<T_y, T_scale, T_shape> loglogistic_cdf(const T_y& y,
                                                     const T_scale& alpha,
                                                     const T_shape& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_scale>;
  using T_beta_ref = ref_type_t<T_shape>;
  using std::pow;
  static constexpr const char* function = "loglogistic_cdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         alpha, "Shape parameter", beta);
  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_nonnegative(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", alpha_val);
  check_positive_finite(function, "Shape parameter", beta_val);

  if (size_zero(y, alpha, beta)) {
    return 1.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  if (sum(promote_scalar<int>(y_val == 0))) {
    return ops_partials.build(0.0);
  }

  const auto& alpha_div_y
      = to_ref_if<!is_constant_all<T_shape>::value>(alpha_val / y_val);
  const auto& alpha_div_y_pow_beta
      = to_ref_if<!is_constant_all<T_y, T_scale, T_shape>::value>(
          pow(alpha_div_y, beta_val));
  const auto& prod_all
      = to_ref_if<!is_constant_all<T_y, T_scale, T_shape>::value>(
          1 / (1 + alpha_div_y_pow_beta));

  T_partials_return cdf = prod(prod_all);

  if (!is_constant_all<T_y, T_scale, T_shape>::value) {
    const auto& prod_all_sq = to_ref_if<!is_constant_all<T_y>::value
                                            + !is_constant_all<T_scale>::value
                                            + !is_constant_all<T_shape>::value
                                        >= 2>(square(prod_all));
    const auto& cdf_div_elt = to_ref_if<!is_constant_all<T_y>::value
                                            + !is_constant_all<T_scale>::value
                                            + !is_constant_all<T_shape>::value
                                        >= 2>(cdf / prod_all);
    if (!is_constant_all<T_y, T_scale>::value) {
      const auto& alpha_div_times_beta = to_ref_if<
          !is_constant_all<T_y>::value + !is_constant_all<T_scale>::value == 2>(
          alpha_div_y_pow_beta * beta_val);
      if (!is_constant_all<T_y>::value) {
        const auto& y_deriv = alpha_div_times_beta / y_val * prod_all_sq;
        partials<0>(ops_partials) = y_deriv * cdf_div_elt;
      }
      if (!is_constant_all<T_scale>::value) {
        const auto& alpha_deriv
            = -alpha_div_times_beta / alpha_val * prod_all_sq;
        partials<1>(ops_partials) = alpha_deriv * cdf_div_elt;
      }
    }
    if (!is_constant_all<T_shape>::value) {
      const auto& beta_deriv
          = -multiply_log(alpha_div_y_pow_beta, alpha_div_y) * prod_all_sq;
      partials<2>(ops_partials) = beta_deriv * cdf_div_elt;
    }
  }

  return ops_partials.build(cdf);
}

}  // namespace math
}  // namespace stan
#endif
