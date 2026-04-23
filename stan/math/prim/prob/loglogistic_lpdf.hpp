#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_LPDF_HPP

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
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the loglogistic density for the specified scalar(s)
 * given the specified scales(s) and shape(s). y, alpha, or beta
 * can each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/scale/shape triple.
 *
 * @tparam T_y type of scalar.
 * @tparam T_scale type of scale parameter.
 * @tparam T_shape type of shape parameter.
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) scale parameter(s)
 * for the loglogistic distribution.
 * @param beta (Sequence of) shape parameter(s) for the
 * loglogistic distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if any of the inputs are not positive and finite.
 */
template <bool propto, typename T_y, typename T_scale, typename T_shape,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale, T_shape>* = nullptr>
inline return_type_t<T_y, T_scale, T_shape> loglogistic_lpdf(
    const T_y& y, const T_scale& alpha, const T_shape& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale, T_shape>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_scale_ref = ref_type_if_not_constant_t<T_scale>;
  using T_shape_ref = ref_type_if_not_constant_t<T_shape>;
  using std::pow;
  static constexpr const char* function = "loglogistic_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Scale parameter",
                         alpha, "Shape parameter", beta);

  T_y_ref y_ref = y;
  T_scale_ref alpha_ref = alpha;
  T_shape_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_positive_finite(function, "Random variable", y_val);
  check_positive_finite(function, "Scale parameter", alpha_val);
  check_positive_finite(function, "Shape parameter", beta_val);

  if (size_zero(y, alpha, beta)) {
    return 0.0;
  }
  if constexpr (!include_summand<propto, T_y, T_scale, T_shape>::value) {
    return 0.0;
  }

  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);

  const auto& inv_alpha
      = to_ref_if<is_any_autodiff_v<T_y, T_scale>>(inv(alpha_val));
  const auto& y_div_alpha
      = to_ref_if<is_autodiff_v<T_shape>>(y_val * inv_alpha);
  const auto& y_div_alpha_pow_beta
      = to_ref_if<is_autodiff_v<T_shape>>(pow(y_div_alpha, beta_val));
  const auto& log1_arg = to_ref_if<is_any_autodiff_v<T_y, T_scale, T_shape>>(
      1 + y_div_alpha_pow_beta);
  const auto& log_y = to_ref_if<is_autodiff_v<T_shape>>(log(y_val));
  const auto& log_alpha
      = to_ref_if<include_summand<propto, T_scale, T_shape>::value>(
          log(alpha_val));
  const auto& beta_minus_one
      = to_ref_if<(include_summand<propto, T_scale, T_shape>::value
                   || is_autodiff_v<T_y>)>(beta_val - 1.0);

  size_t N = max_size(y, alpha, beta);
  size_t N_alpha_beta = max_size(alpha, beta);

  T_partials_return logp = sum(beta_minus_one * log_y - 2.0 * log(log1_arg));

  if constexpr (include_summand<propto, T_scale, T_shape>::value) {
    logp += sum(N * (log(beta_val) - log_alpha - beta_minus_one * log_alpha)
                / N_alpha_beta);
  }

  if constexpr (is_any_autodiff_v<T_y, T_scale, T_shape>) {
    const auto& two_inv_log1_arg = to_ref_if<
        is_autodiff_v<
            T_y> + is_autodiff_v<T_scale> + is_autodiff_v<T_shape> >= 2>(
        2.0 * inv(log1_arg));
    if constexpr (is_any_autodiff_v<T_y, T_scale>) {
      const auto& y_pow_beta
          = to_ref_if<is_any_autodiff_v<T_y, T_scale>>(pow(y_val, beta_val));
      const auto& inv_alpha_pow_beta
          = to_ref_if < is_autodiff_v<
                T_y> && is_autodiff_v<T_scale> > (pow(inv_alpha, beta_val));

      if constexpr (is_autodiff_v<T_y>) {
        const auto& inv_y = inv(y_val);
        const auto& y_deriv = beta_minus_one * inv_y
                              - two_inv_log1_arg
                                    * (beta_val * inv_alpha_pow_beta)
                                    * y_pow_beta * inv_y;
        partials<0>(ops_partials) = y_deriv;
      }
      if constexpr (is_autodiff_v<T_scale>) {
        const auto& alpha_deriv = -beta_val * inv_alpha
                                  - two_inv_log1_arg * y_pow_beta * (-beta_val)
                                        * inv_alpha_pow_beta * inv_alpha;
        partials<1>(ops_partials) = alpha_deriv;
      }
    }
    if constexpr (is_autodiff_v<T_shape>) {
      const auto& beta_deriv
          = (1.0 * inv(beta_val)) + log_y - log_alpha
            - two_inv_log1_arg * y_div_alpha_pow_beta * log(y_div_alpha);
      partials<2>(ops_partials) = beta_deriv;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_scale, T_shape> loglogistic_lpdf(
    const T_y& y, const T_scale& alpha, const T_shape& beta) {
  return loglogistic_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
