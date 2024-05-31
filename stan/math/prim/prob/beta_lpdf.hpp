#ifndef STAN_MATH_PRIM_PROB_BETA_LPDF_HPP
#define STAN_MATH_PRIM_PROB_BETA_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the beta density for the specified scalar(s) given the specified
 * sample stan::math::size(s). y, alpha, or beta can each either be scalar or a
 * vector. Any vector inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/alpha/beta triple.
 *
 * Prior sample sizes, alpha and beta, must be greater than 0.
 *
 * @tparam T_y type of scalar outcome
 * @tparam T_scale_succ type of prior scale for successes
 * @tparam T_scale_fail type of prior scale for failures
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) prior sample stan::math::size(s).
 * @param beta (Sequence of) prior sample stan::math::size(s).
 * @return The log of the product of densities.
 */
template <bool propto, typename T_y, typename T_scale_succ,
          typename T_scale_fail,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_scale_succ, T_scale_fail>* = nullptr>
return_type_t<T_y, T_scale_succ, T_scale_fail> beta_lpdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  using T_partials_return = partials_return_t<T_y, T_scale_succ, T_scale_fail>;
  using T_y_ref = ref_type_if_not_constant_t<T_y>;
  using T_alpha_ref = ref_type_if_not_constant_t<T_scale_succ>;
  using T_beta_ref = ref_type_if_not_constant_t<T_scale_fail>;
  static constexpr const char* function = "beta_lpdf";
  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta);
  if (size_zero(y, alpha, beta)) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) alpha_val = to_ref(as_value_column_array_or_scalar(alpha_ref));
  decltype(auto) beta_val = to_ref(as_value_column_array_or_scalar(beta_ref));

  check_positive_finite(function, "First shape parameter", alpha_val);
  check_positive_finite(function, "Second shape parameter", beta_val);
  check_bounded(function, "Random variable", value_of(y_val), 0, 1);
  if (!include_summand<propto, T_y, T_scale_succ, T_scale_fail>::value) {
    return 0;
  }

  const auto& log_y = to_ref(log(y_val));
  const auto& log1m_y = to_ref(log1m(y_val));

  size_t N = max_size(y, alpha, beta);
  T_partials_return logp(0);
  if (include_summand<propto, T_scale_succ>::value) {
    logp -= sum(lgamma(alpha_val)) * N / max_size(alpha);
  }
  if (include_summand<propto, T_scale_fail>::value) {
    logp -= sum(lgamma(beta_val)) * N / max_size(beta);
  }
  if (include_summand<propto, T_y, T_scale_succ>::value) {
    logp += sum((alpha_val - 1.0) * log_y) * N / max_size(y, alpha);
  }
  if (include_summand<propto, T_y, T_scale_fail>::value) {
    logp += sum((beta_val - 1.0) * log1m_y) * N / max_size(y, beta);
  }

  auto ops_partials = make_partials_propagator(y_ref, alpha_ref, beta_ref);
  if (!is_constant_all<T_y>::value) {
    edge<0>(ops_partials).partials_
        = (alpha_val - 1) / y_val + (beta_val - 1) / (y_val - 1);
  }

  if (include_summand<propto, T_scale_succ, T_scale_fail>::value) {
    const auto& alpha_beta
        = to_ref_if<!is_constant_all<T_scale_succ, T_scale_fail>::value>(
            alpha_val + beta_val);
    logp += sum(lgamma(alpha_beta)) * N / max_size(alpha, beta);
    if (!is_constant_all<T_scale_succ, T_scale_fail>::value) {
      const auto& digamma_alpha_beta
          = to_ref_if < !is_constant_all<T_scale_succ>::value
            && !is_constant_all<T_scale_fail>::value > (digamma(alpha_beta));
      if (!is_constant_all<T_scale_succ>::value) {
        edge<1>(ops_partials).partials_
            = log_y + digamma_alpha_beta - digamma(alpha_val);
      }
      if (!is_constant_all<T_scale_fail>::value) {
        edge<2>(ops_partials).partials_
            = log1m_y + digamma_alpha_beta - digamma(beta_val);
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline return_type_t<T_y, T_scale_succ, T_scale_fail> beta_lpdf(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  return beta_lpdf<false>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
