#ifndef STAN_MATH_PRIM_PROB_WEIBULL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WEIBULL_LPDF_HPP

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
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the Weibull log probability density for the given
 * location and scale. Given containers of matching sizes, returns the
 * log sum of probability densities.
 *
 * @tparam T_y type of real parameter
 * @tparam T_shape type of shape parameter
 * @tparam T_scale type of scale parameter
 * @param y real parameter
 * @param alpha shape parameter
 * @param sigma scale parameter
 * @return log probability density or log sum of probability densities
 * @throw std::domain_error if y is negative, alpha or sigma are nonpositive
 */
template <bool propto, typename T_y, typename T_shape, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_scale>* = nullptr>
return_type_t<T_y, T_shape, T_scale> weibull_lpdf(const T_y& y,
                                                  const T_shape& alpha,
                                                  const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  using std::pow;
  static const char* function = "weibull_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", sigma);

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_finite(function, "Random variable", y_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Scale parameter", sigma_val);

  if (size_zero(y, alpha, sigma)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_shape, T_scale>::value) {
    return 0;
  }

  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, sigma_ref);

  if (sum(promote_scalar<int>(y_val < 0))) {
    return LOG_ZERO;
  }

  const auto& log_y
      = to_ref_if<include_summand<propto, T_y, T_shape>::value>(log(y_val));
  const auto& log_sigma
      = to_ref_if<include_summand<propto, T_shape, T_scale>::value>(
          log(sigma_val));
  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_scale>::value>(inv(sigma_val));
  const auto& y_div_sigma_pow_alpha
      = to_ref_if<!is_constant_all<T_y, T_shape, T_scale>::value>(
          pow(y_val * inv_sigma, alpha_val));

  size_t N = max_size(y, alpha, sigma);
  T_partials_return logp = -sum(y_div_sigma_pow_alpha);
  if (include_summand<propto, T_shape>::value) {
    logp += sum(log(alpha_val)) * N / math::size(alpha);
  }
  if (include_summand<propto, T_y, T_shape>::value) {
    logp += sum((alpha_val - 1.0) * log_y) * N / max_size(alpha, y);
  }
  if (include_summand<propto, T_shape, T_scale>::value) {
    logp -= sum(alpha_val * log_sigma) * N / max_size(alpha, sigma);
  }

  if (!is_constant_all<T_y>::value) {
    edge<0>(ops_partials).partials_
        = (alpha_val * (1 - y_div_sigma_pow_alpha) - 1.0) / y_val;
  }
  if (!is_constant_all<T_shape>::value) {
    edge<1>(ops_partials).partials_
        = inv(alpha_val) + (1.0 - y_div_sigma_pow_alpha) * (log_y - log_sigma);
  }
  if (!is_constant_all<T_scale>::value) {
    edge<2>(ops_partials).partials_
        = alpha_val * inv_sigma * (y_div_sigma_pow_alpha - 1.0);
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_shape, typename T_scale>
inline return_type_t<T_y, T_shape, T_scale> weibull_lpdf(const T_y& y,
                                                         const T_shape& alpha,
                                                         const T_scale& sigma) {
  return weibull_lpdf<false>(y, alpha, sigma);
}

}  // namespace math
}  // namespace stan
#endif
