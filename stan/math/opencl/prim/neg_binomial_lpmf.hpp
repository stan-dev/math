#ifndef STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_LPMF_HPP
#define STAN_MATH_OPENCL_PRIM_NEG_BINOMIAL_LPMF_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the negative binomial density for the specified scalars given
 * the specified mean(s) and deviation(s). n, alpha, or beta can
 * each be either a scalar or a vector matrix_cl. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_n_cl type of scalar
 * @tparam T_shape_cl type of location parameter
 * @tparam T_inv_scale_cl type of precision parameter
 * @param n (Sequence of) scalar(s).
 * @param alpha (Sequence of) location parameter(s)
 * @param beta (Sequence of) precision parameters
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, typename T_n_cl, typename T_shape_cl,
          typename T_inv_scale_cl,
          require_all_prim_or_rev_kernel_expression_t<
              T_n_cl, T_shape_cl, T_inv_scale_cl>* = nullptr,
          require_any_not_stan_scalar_t<T_n_cl, T_shape_cl,
                                        T_inv_scale_cl>* = nullptr>
inline return_type_t<T_n_cl, T_shape_cl, T_inv_scale_cl> neg_binomial_lpmf(
    const T_n_cl& n, const T_shape_cl& alpha, const T_inv_scale_cl& beta) {
  static const char* function = "neg_binomial_lpmf(OpenCL)";
  using T_partials_return
      = partials_return_t<T_n_cl, T_shape_cl, T_inv_scale_cl>;
  using std::isfinite;
  using std::isnan;

  check_consistent_sizes(function, "Failures variable", n, "Shape parameter",
                         alpha, "Inverse scale parameter", beta);
  const size_t N = max_size(n, alpha, beta);
  if (N == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_n_cl, T_shape_cl, T_inv_scale_cl>::value) {
    return 0.0;
  }

  const auto& alpha_col = as_column_vector_or_scalar(alpha);
  const auto& beta_col = as_column_vector_or_scalar(beta);

  const auto& alpha_val = value_of(alpha_col);
  const auto& beta_val = value_of(beta_col);

  auto check_n_nonnegative
      = check_cl(function, "Failures variable", n, "nonnegative");
  auto n_nonnegative = n >= 0;
  auto check_alpha_positive_finite
      = check_cl(function, "Shape parameter", alpha_val, "positive finite");
  auto alpha_positive_finite = 0 < alpha_val && isfinite(alpha_val);
  auto check_beta_positive_finite = check_cl(
      function, "Inverse scale parameter", beta_val, "positive finite");
  auto beta_positive_finite = 0 < beta_val && isfinite(beta_val);

  auto digamma_alpha = digamma(alpha_val);
  auto log1p_inv_beta = log1p(elt_divide(1.0, beta_val));
  auto log1p_beta = log1p(beta_val);
  auto lambda_m_alpha_over_1p_beta
      = elt_divide(alpha_val, beta_val) - elt_divide(alpha_val, 1.0 + beta_val);

  auto logp1
      = -elt_multiply(alpha_val, log1p_inv_beta) - elt_multiply(n, log1p_beta);
  auto logp_expr
      = colwise_sum(static_select<include_summand<propto, T_shape_cl>::value>(
          logp1
              + binomial_coefficient_log(n + alpha_val - 1.0, alpha_val - 1.0),
          logp1));

  auto alpha_deriv = digamma(alpha_val + n) - digamma_alpha - log1p_inv_beta;
  auto beta_deriv = lambda_m_alpha_over_1p_beta - elt_divide(n, beta_val + 1.0);

  matrix_cl<double> logp_cl;
  matrix_cl<double> alpha_deriv_cl;
  matrix_cl<double> beta_deriv_cl;

  results(check_n_nonnegative, check_alpha_positive_finite,
          check_beta_positive_finite, logp_cl, alpha_deriv_cl, beta_deriv_cl)
      = expressions(n_nonnegative, alpha_positive_finite, beta_positive_finite,
                    logp_expr,
                    calc_if<!is_constant<T_shape_cl>::value>(alpha_deriv),
                    calc_if<!is_constant<T_inv_scale_cl>::value>(beta_deriv));

  T_partials_return logp = sum(from_matrix_cl(logp_cl));

  operands_and_partials<decltype(alpha_col), decltype(beta_col)> ops_partials(
      alpha_col, beta_col);

  if (!is_constant<T_shape_cl>::value) {
    ops_partials.edge1_.partials_ = std::move(alpha_deriv_cl);
  }
  if (!is_constant<T_inv_scale_cl>::value) {
    ops_partials.edge2_.partials_ = std::move(beta_deriv_cl);
  }
  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
