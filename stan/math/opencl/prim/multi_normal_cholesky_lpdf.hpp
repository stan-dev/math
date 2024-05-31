#ifndef STAN_MATH_OPENCL_PRIM_MULTI_NORMAL_CHOLESKY_LPDF_HPP
#define STAN_MATH_OPENCL_PRIM_MULTI_NORMAL_CHOLESKY_LPDF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/dot_self.hpp>
#include <stan/math/opencl/prim/mdivide_left_tri_low.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>
#include <stan/math/prim/fun/elt_multiply.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <stan/math/prim/err/constraint_tolerance.hpp>

namespace stan {
namespace math {

/** \ingroup opencl
 * The log of the multivariate normal density for the given y, mu, and
 * a Cholesky factor L of the variance matrix.
 * Sigma = LL', a square, semi-positive definite matrix.
 *
 * Analytic expressions taken from
 * http://qwone.com/~jason/writing/multivariateNormal.pdf
 * written by Jason D. M. Rennie.
 *
 * @param y A scalar vector
 * @param mu The mean vector of the multivariate normal distribution.
 * @param L The Cholesky decomposition of a variance matrix
 * of the multivariate normal distribution
 * @return The log of the multivariate normal density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_loc Type of location.
 * @tparam T_covar Type of scale.
 */
template <bool propto, typename T_y_cl, typename T_loc_cl, typename T_covar_cl,
          require_all_prim_or_rev_kernel_expression_t<T_y_cl, T_loc_cl,
                                                      T_covar_cl>* = nullptr>
inline return_type_t<T_y_cl, T_loc_cl, T_covar_cl> multi_normal_cholesky_lpdf(
    const T_y_cl& y, const T_loc_cl& mu, const T_covar_cl& L) {
  static constexpr const char* function = "multi_normal_cholesky_lpdf(OpenCL)";

  check_consistent_sizes(function, "y", y, "mu", mu);
  check_square(function, "covariance parameter", L);
  check_size_match(function, "Size of random variable", y.rows(),
                   "rows of covariance parameter", L.rows());

  if (max_size(y, mu, L) == 0) {
    return 0.0;
  }
  if (!include_summand<propto, T_y_cl, T_loc_cl, T_covar_cl>::value) {
    return 0.0;
  }

  const auto& y_val = value_of(y);
  const auto& mu_val = value_of(mu);
  const auto& L_val_eval = eval(value_of(L));

  int L_size = L_val_eval.rows();
  int N_cases = std::max(y_val.cols(), mu_val.cols());

  double logp = 0;
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * L_size * N_cases;
  }

  matrix_cl<double> L_lower(L_val_eval.buffer(), L_val_eval.rows(),
                            L_val_eval.cols(), matrix_cl_view::Lower);
  matrix_cl<double> inv_L = mdivide_left_tri_low(L_lower);

  auto check_mu_finite
      = check_cl(function, "Location parameter", mu_val, "finite");
  auto mu_finite = isfinite(mu_val);
  auto check_y_not_nan
      = check_cl(function, "Random variable", y_val, "not nan");
  auto y_not_nan = !isnan(y_val);

  auto sum_log_diag_inv_L = colwise_sum(log(diagonal(inv_L)));
  auto y_mu_diff
      = rowwise_optional_broadcast(y_val) - rowwise_optional_broadcast(mu_val);

  matrix_cl<double> y_mu_diff_cl;
  matrix_cl<double> sum_log_diag_inv_L_cl;

  if (y_val.cols() == 1 && mu_val.cols() == 1) {
    results(check_mu_finite, check_y_not_nan, y_mu_diff_cl,
            sum_log_diag_inv_L_cl)
        = expressions(mu_finite, y_not_nan, y_mu_diff,
                      calc_if<include_summand<propto, T_covar_cl>::value>(
                          sum_log_diag_inv_L));
  } else if (y_val.cols() == 1) {
    results(check_y_not_nan, sum_log_diag_inv_L_cl) = expressions(
        y_not_nan, calc_if<include_summand<propto, T_covar_cl>::value>(
                       sum_log_diag_inv_L));
    results(check_mu_finite, y_mu_diff_cl) = expressions(mu_finite, y_mu_diff);
  } else if (mu_val.cols() == 1) {
    results(check_mu_finite, sum_log_diag_inv_L_cl) = expressions(
        mu_finite, calc_if<include_summand<propto, T_covar_cl>::value>(
                       sum_log_diag_inv_L));
    results(check_y_not_nan, y_mu_diff_cl) = expressions(y_not_nan, y_mu_diff);
  } else {
    sum_log_diag_inv_L_cl = calc_if<include_summand<propto, T_covar_cl>::value>(
        sum_log_diag_inv_L);
    results(check_mu_finite, check_y_not_nan, y_mu_diff_cl)
        = expressions(mu_finite, y_not_nan, y_mu_diff);
  }

  if (include_summand<propto, T_covar_cl>::value) {
    logp += sum(from_matrix_cl(sum_log_diag_inv_L_cl)) * N_cases;
  }

  matrix_cl<double> half = transpose(inv_L * y_mu_diff_cl);
  matrix_cl<double> scaled_diff = transpose(half * inv_L);
  logp -= 0.5 * dot_self(half);

  auto ops_partials = make_partials_propagator(y, mu, L);

  if (!is_constant_all<T_y_cl>::value) {
    if (y_val.cols() == 1) {
      partials<0>(ops_partials) = -rowwise_sum(scaled_diff);
    } else {
      partials<0>(ops_partials) = -scaled_diff;
    }
  }
  if (!is_constant_all<T_loc_cl>::value) {
    if (mu_val.cols() == 1) {
      partials<1>(ops_partials) = rowwise_sum(scaled_diff);
    } else {
      partials<1>(ops_partials) = scaled_diff;
    }
  }
  if (!is_constant_all<T_covar_cl>::value) {
    partials<2>(ops_partials) = scaled_diff * half - N_cases * transpose(inv_L);
  }

  return ops_partials.build(logp);
}

}  // namespace math
}  // namespace stan
#endif
#endif
