#ifndef STAN_MATH_OPENCL_PRIM_MULTI_NORMAL_CHOLESKY_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_MULTI_NORMAL_CHOLESKY_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/mdivide_left_tri_low.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * Returns the (natural) log probability of the specified array
 * of integers given the vector of continuous locations and
 * specified cutpoints in an ordered logistic model.
 *
 * <p>Typically the continuous location
 * will be the dot product of a vector of regression coefficients
 * and a vector of predictors for the outcome
 *
  \f[
    \frac{\partial }{\partial \lambda} =
    \begin{cases}\\
    -\mathrm{logit}^{-1}(\lambda - c_1) & \mbox{if } k = 1,\\
    -(((1-e^{c_{k-1}-c_{k-2}})^{-1} - \mathrm{logit}^{-1}(c_{k-2}-\lambda)) +
    ((1-e^{c_{k-2}-c_{k-1}})^{-1} - \mathrm{logit}^{-1}(c_{k-1}-\lambda)))
    & \mathrm{if } 1 < k < K, \mathrm{and}\\
    \mathrm{logit}^{-1}(c_{K-2}-\lambda) & \mathrm{if } k = K.
    \end{cases}
  \f]

  \f[
    \frac{\partial }{\partial \lambda} =
    \begin{cases}
    -\mathrm{logit}^{-1}(\lambda - c_1) & \text{if } k = 1,\\
    -(((1-e^{c_{k-1}-c_{k-2}})^{-1} - \mathrm{logit}^{-1}(c_{k-2}-\lambda)) +
    ((1-e^{c_{k-2}-c_{k-1}})^{-1} - \mathrm{logit}^{-1}(c_{k-1}-\lambda)))
    & \text{if } 1 < k < K, \text{ and}\\
    \mathrm{logit}^{-1}(c_{K-2}-\lambda) & \text{if } k = K.
    \end{cases}
  \f]
 *
 * @tparam propto True if calculating up to a proportion.
 * @tparam T_y Y variable type (integer or array of integers).
 * @tparam T_loc Location type.
 * @tparam T_cut Cut-point type.
 * @param y Array of integers
 * @param lambda Vector of continuous location variables.
 * @param c Positive increasing vector of cutpoints.
 * @return Log probability of outcome given location and
 * cutpoints.
 * @throw std::domain_error If the outcome is not between 1 and
 * the number of cutpoints plus 2; if the cutpoint vector is
 * empty; if the cutpoint vector contains a non-positive,
 * non-finite value; or if the cutpoint vector is not sorted in
 * ascending order.
 * @throw std::invalid_argument If y and lambda are different
 * lengths.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_cuts_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                      T_cuts_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_cut_cl> ordered_logistic_lpmf(
    const T_y_cl& y, const T_loc_cl& lambda, const T_cuts_cl& cuts) {
  static const char* function = "ordered_logistic_lpmf(OpenCL)";

  check_nonzero_size(function, "Cut-points", c);
  check_consistent_sizes(function, "Integers", y, "Locations", lambda);

  int N_instances = max_size(y, lambda);
  int N_classes = cuts.rows() + 1;
  int N_cut_sets = cuts.cols();

  if (N_cut_sets > 1) {
    check_size_match(function, "Length of location variables ", N_instances,
                     "Number of cutpoint vectors ", N_cut_sets);
  }

  if (max_size(y, lambda) == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_cuts_cl>::value) {
    return 0.0;
  }

  const auto& y_val = value_of(y);
  const auto& lambda_val = value_of(lambda);
  const auto& cuts_val = value_of(cuts);

  if (N_classes >= 2) {
    auto cuts_head
        = block_zero_based(cuts_val, 0, 0, size(cuts) - 1, N_cut_sets);
    auto cuts_tail
        = block_zero_based(cuts_val, 1, 0, size(cuts) - 1, N_cut_sets);
    check_cl(function, "Cuts", cuts_head, "ordered and finite")
        = cuts_head < cuts_tail && isfinite(cuts_head) && isfinite(cuts_tail);
  } else {
    check_cl(function, "Cuts", cuts_val, "finite") = isfinite(cuts_val);
  }

  auto cuts_y1 = select(y_val != N_classes,
                        indexing(cuts_val, y_val - 1, col_index()), INFINITY);
  auto cuts_y2 = select(y_val != 1, indexing(cuts_val, y_val - 2, col_index()),
                        -INFINITY);
  auto cut2 = lambda_val - cuts_y2;
  auto cut1 = lambda_val - cuts_y1;
  auto m_log_1p_exp_cut1 = -log1p_exp(cut1);
  auto m_log_1p_exp_m_cut2 = -log1p_exp(-cut2);
  auto log1m_exp_cuts_diff = log1m_exp(cut1 - cut2);
  auto logp_tmp = colwise_sum(
      y_val == 1, m_log_1p_exp_cut1,
      select(y_val == N_classes, m_log_1p_exp_m_cut2,
             m_log_1p_exp_m_cut2 + log1m_exp_cuts_diff + m_log_1p_exp_cut1));

  auto exp_m_cut1 = exp(-cut1);
  auto exp_m_cut2 = exp(-cut2);
  auto exp_cuts_diff = exp(cuts_y2 - cuts_y1);
  auto d1
      = select(cut2 > 0.0, exp_m_cut2 / (1.0 + exp_m_cut2), 1.0 / (1.0 + exp(cut2)))
        - exp_cuts_diff / (exp_cuts_diff - 1.0);
  auto d2
      = 1.0 / (1.0 - exp_cuts_diff)
        - select(cut1 > 0.0, exp_m_cut1 / (1.0 + exp_m_cut1),
                            1.0 / (1.0 + exp(cut1)));


  operands_and_partials<T_y_cl, T_loc_cl, T_cuts_cl> ops_partials(y, lambda,
                                                                   cuts);



  int L_size = L_val_eval.rows();
  int N_cases = std::max(y_val.cols(), lambda_val.cols());

  double logp = 0;
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * L_size * N_cases;
  }

  matrix_cl<double> L_lower(L_val_eval.buffer(), L_val_eval.rows(),
                            L_val_eval.cols(), matrix_cl_view::Lower);
  matrix_cl<double> inv_L = mdivide_left_tri_low(L_lower);

  auto check_lambda_finite
      = check_cl(function, "Location parameter", lambda_val, "finite");
  auto lambda_finite = isfinite(lambda_val);
  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not nan");
  auto y_not_nan = !isnan(y_val);

  auto sum_log_diag_inv_L = colwise_sum(log(diagonal(inv_L)));
  auto y_lambda_diff = rowwise_optional_broadcast(y_val)
                       - rowwise_optional_broadcast(lambda_val);

  matrix_cl<double> y_lambda_diff_cl;
  matrix_cl<double> sum_log_diag_inv_L_cl;

  if (y_val.cols() == 1 && lambda_val.cols() == 1) {
    results(check_lambda_finite, check_y_not_nan, y_lambda_diff_cl,
            sum_log_diag_inv_L_cl)
        = expressions(lambda_finite, y_not_nan, y_lambda_diff,
                      calc_if<include_summand<propto, T_covar_cl>::value>(
                          sum_log_diag_inv_L));
  } else if (y_val.cols() == 1) {
    results(check_y_not_nan, sum_log_diag_inv_L_cl) = expressions(
        y_not_nan, calc_if<include_summand<propto, T_covar_cl>::value>(
                       sum_log_diag_inv_L));
    results(check_lambda_finite, y_lambda_diff_cl)
        = expressions(lambda_finite, y_lambda_diff);
  } else if (lambda_val.cols() == 1) {
    results(check_lambda_finite, sum_log_diag_inv_L_cl) = expressions(
        lambda_finite, calc_if<include_summand<propto, T_covar_cl>::value>(
                           sum_log_diag_inv_L));
    results(check_y_not_nan, y_lambda_diff_cl)
        = expressions(y_not_nan, y_lambda_diff);
  } else {
    sum_log_diag_inv_L_cl = calc_if<include_summand<propto, T_covar_cl>::value>(
        sum_log_diag_inv_L);
    results(check_lambda_finite, check_y_not_nan, y_lambda_diff_cl)
        = expressions(lambda_finite, y_not_nan, y_lambda_diff);
  }

  if (include_summand<propto, T_covar_cl>::value) {
    logp += sum(from_matrix_cl(sum_log_diag_inv_L_cl)) * N_cases;
  }

  matrix_cl<double> half = transpose(inv_L * y_lambda_diff_cl);
  matrix_cl<double> scaled_diff = transpose(half * inv_L);
  logp -= 0.5 * dot_self(half);

  operands_and_partials<T_y_cl, T_loc_cl, T_covar_cl> ops_partials(y, lambda,
                                                                   L);

  if (!is_constant_all<T_y_cl>::value) {
    if (y_val.cols() == 1) {
      forward_as<matrix_cl<double>>(ops_partials.edge1_.partials_)
          = -rowwise_sum(scaled_diff);
    } else {
      ops_partials.edge1_.partials_ = -scaled_diff;
    }
  }
  if (!is_constant_all<T_loc_cl>::value) {
    if (lambda_val.cols() == 1) {
      forward_as<matrix_cl<double>>(ops_partials.edge2_.partials_)
          = rowwise_sum(scaled_diff);
    } else {
      ops_partials.edge2_.partials_ = scaled_diff;
    }
  }
  if (!is_constant_all<T_covar_cl>::value) {
    ops_partials.edge3_.partials_
        = scaled_diff * half - N_cases * transpose(inv_L);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
