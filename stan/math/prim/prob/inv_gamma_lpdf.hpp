#ifndef STAN_MATH_PRIM_PROB_INV_GAMMA_LPDF_HPP
#define STAN_MATH_PRIM_PROB_INV_GAMMA_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
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

/** \ingroup prob_dists
 * The log of an inverse gamma density for y with the specified
 * shape and scale parameters.
 * Shape and scale parameters must be greater than 0.
 * y must be greater than 0.
 *
 * @param y A scalar variable.
 * @param alpha Shape parameter.
 * @param beta Scale parameter.
 * @throw std::domain_error if alpha is not greater than 0.
 * @throw std::domain_error if beta is not greater than 0.
 * @throw std::domain_error if y is not greater than 0.
 * @tparam T_y Type of scalar.
 * @tparam T_shape Type of shape.
 * @tparam T_scale Type of scale.
 */
template <bool propto, typename T_y, typename T_shape, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_shape, T_scale>* = nullptr>
return_type_t<T_y, T_shape, T_scale> inv_gamma_lpdf(const T_y& y,
                                                    const T_shape& alpha,
                                                    const T_scale& beta) {
  using T_partials_return = partials_return_t<T_y, T_shape, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_shape>::value, T_shape>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "inv_gamma_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Shape parameter",
                         alpha, "Scale parameter", beta);

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_not_nan(function, "Random variable", y_val);
  check_positive_finite(function, "Shape parameter", alpha_val);
  check_positive_finite(function, "Scale parameter", beta_val);

  if (size_zero(y, alpha, beta)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_shape, T_scale>::value) {
    return 0;
  }
  if (sum(promote_scalar<int>(y_val <= 0))) {
    return LOG_ZERO;
  }

  T_partials_return logp(0);
  operands_and_partials<T_y_ref, T_alpha_ref, T_beta_ref> ops_partials(
      y_ref, alpha_ref, beta_ref);

  const auto& log_y
      = to_ref_if<include_summand<propto, T_y, T_shape>::value>(log(y_val));

  size_t N = max_size(y, alpha, beta);
  if (include_summand<propto, T_shape>::value) {
    logp -= sum(lgamma(alpha_val)) * N / size(alpha);
  }
  if (include_summand<propto, T_shape, T_scale>::value) {
    const auto& log_beta
        = to_ref_if<!is_constant_all<T_shape>::value>(log(beta_val));
    logp += sum(alpha_val * log_beta) * N / max_size(alpha, beta);
    if (!is_constant_all<T_shape>::value) {
      ops_partials.edge2_.partials_ = log_beta - digamma(alpha_val) - log_y;
    }
  }
  if (include_summand<propto, T_y, T_shape>::value) {
    logp -= sum((alpha_val + 1.0) * log_y) * N / max_size(y, alpha);
  }
  if (include_summand<propto, T_y, T_scale>::value) {
    const auto& inv_y
        = to_ref_if<(!is_constant_all<T_y>::value
                     || !is_constant_all<T_scale>::value)>(inv(y_val));
    logp -= sum(beta_val * inv_y) * N / max_size(y, beta);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_
          = (beta_val * inv_y - alpha_val - 1) * inv_y;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_ = alpha_val / beta_val - inv_y;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_shape, typename T_scale>
inline return_type_t<T_y, T_shape, T_scale> inv_gamma_lpdf(
    const T_y& y, const T_shape& alpha, const T_scale& beta) {
  return inv_gamma_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
