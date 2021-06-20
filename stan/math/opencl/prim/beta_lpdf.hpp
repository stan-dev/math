#ifndef STAN_MATH_OPENCL_PRIM_BETA_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_BETA_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/size.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the beta density for the specified scalar(s) given the specified
 * sample stan::math::size(s). y, alpha, or beta can each either be scalar or a
 * vector on OpenCL device. Any vector inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/alpha/beta triple.
 *
 * Prior sample sizes, alpha and beta, must be greater than 0.
 *
 * @tparam T_y_cl type of scalar outcome
 * @tparam T_scale_succ_cl type of prior scale for successes
 * @tparam T_scale_fail_cl type of prior scale for failures
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) prior sample stan::math::size(s).
 * @param beta (Sequence of) prior sample stan::math::size(s).
 * @return The log of the product of densities.
 */
template <bool propto, typename T_y_cl, typename T_scale_succ_cl,
          typename T_scale_fail_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_y_cl, T_scale_succ_cl, T_scale_fail_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_y_cl, T_scale_succ_cl,
                                        T_scale_fail_cl>* = nullptr>
return_type_t<T_y_cl, T_scale_succ_cl, T_scale_fail_cl> beta_lpdf(
    const T_y_cl& y, const T_scale_succ_cl& alpha,
    const T_scale_fail_cl& beta) {
  using std::isfinite;
  static const char* function = "beta_lpdf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_y_cl, T_scale_succ_cl, T_scale_fail_cl>;

  check_consistent_sizes(function, "Random variable", y,
                         "First shape parameter", alpha,
                         "Second shape parameter", beta);
  const size_t N = max_size(y, alpha, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_scale_succ_cl,
                       T_scale_fail_cl>::value) {
    return 0.0;
  }

  const auto& y_col = as_column_vector_or_scalar(y);
  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& y_val = value_of(y_col);
  const auto& alpha_val = value_of(alpha_col);
  const auto& beta_val = value_of(beta_col);

  operands_and_partials<decltype(y_col), decltype(alpha_col),
                        decltype(beta_col)>
      ops_partials(y_col, alpha_col, beta_col);

  auto check_alpha_pos_finite = check_cl(function, "First shape parameter",
                                         alpha_val, "positive finite");
  auto alpha_pos_finite = alpha_val > 0 && isfinite(alpha_val);
  auto check_beta_pos_finite = check_cl(function, "Second shape parameter",
                                        beta_val, "positive finite");
  auto beta_pos_finite = beta_val > 0 && isfinite(beta_val);
  auto check_y_bounded
      = check_cl(function, "Random variable", y_val, "in the interval [0, 1]");
  auto y_bounded = 0 <= y_val && y_val <= 1;

  auto log_y_expr = log(y_val);
  auto log1m_y_expr = log1p(-y_val);
  auto alpha_beta_expr = alpha_val + beta_val;

  auto zero_expr
      = as_operation_cl(0);  // simplifiy the kernel by only using one zero
  auto logp_expr = colwise_sum(
      static_select<include_summand<propto, T_scale_succ_cl>::value>(
          -lgamma(alpha_val), zero_expr)
      + static_select<include_summand<propto, T_scale_fail_cl>::value>(
            -lgamma(beta_val), zero_expr)
      + static_select<include_summand<propto, T_y_cl, T_scale_succ_cl>::value>(
            elt_multiply((alpha_val - 1.0), log_y_expr), zero_expr)
      + static_select<include_summand<propto, T_y_cl, T_scale_fail_cl>::value>(
            elt_multiply((beta_val - 1.0), log1m_y_expr), zero_expr)
      + static_select<
            include_summand<propto, T_scale_succ_cl, T_scale_fail_cl>::value>(
            lgamma(alpha_beta_expr), zero_expr));

  auto y_deriv_expr = calc_if<!is_constant<T_y_cl>::value>(
      elt_divide((alpha_val - 1), y_val)
      + elt_divide((beta_val - 1), (y_val - 1)));
  auto digamma_alpha_beta_expr = digamma(alpha_beta_expr);
  auto alpha_deriv_expr = calc_if<!is_constant<T_scale_succ_cl>::value>(
      log_y_expr + digamma_alpha_beta_expr - digamma(alpha_val));
  auto beta_deriv_expr = calc_if<!is_constant<T_scale_fail_cl>::value>(
      log1m_y_expr + digamma_alpha_beta_expr - digamma(beta_val));

  matrix_cl<double> logp_cl;
  matrix_cl<double> y_deriv_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_alpha_pos_finite, check_beta_pos_finite, check_y_bounded,
          logp_cl, y_deriv_cl, alpha_deriv_cl, beta_deriv_cl)
      = expressions(alpha_pos_finite, beta_pos_finite, y_bounded, logp_expr,
                    y_deriv_expr, alpha_deriv_expr, beta_deriv_expr);

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  if (!is_constant<T_y_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(y_deriv_cl);
  }
  if (!is_constant<T_scale_succ_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_scale_fail_cl>::value) {
    ops_partials.edge3_.partials_ = std::move(beta_deriv_cl);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
