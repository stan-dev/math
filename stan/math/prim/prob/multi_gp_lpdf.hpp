#ifndef STAN_MATH_PRIM_PROB_MULTI_GP_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_GP_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/trace_gen_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of a multivariate Gaussian Process for the given y, Sigma, and
 * w.  y is a dxN matrix, where each column is a different observation and each
 * row is a different output dimension.  The Gaussian Process is assumed to
 * have a scaled kernel matrix with a different scale for each output dimension.
 * This distribution is equivalent to:
 *    for (i in 1:d) row(y, i) ~ multi_normal(0, (1/w[i])*Sigma).
 *
 * @tparam T_y type of scalar
 * @tparam T_covar type of kernel
 * @tparam T_w type of weight
 * @param y A dxN matrix
 * @param Sigma The NxN kernel matrix
 * @param w A d-dimensional vector of positive inverse scale parameters for each
 * output.
 * @return The log of the multivariate GP density.
 * @throw std::domain_error if Sigma is not square, not symmetric,
 * or not semi-positive definite.
 */
template <bool propto, typename T_y, typename T_covar, typename T_w,
          require_all_matrix_t<T_y, T_covar>* = nullptr,
          require_col_vector_t<T_w>* = nullptr>
return_type_t<T_y, T_covar, T_w> multi_gp_lpdf(const T_y& y,
                                               const T_covar& Sigma,
                                               const T_w& w) {
  using T_lp = return_type_t<T_y, T_covar, T_w>;
  static const char* function = "multi_gp_lpdf";
  check_size_match(function, "Size of random variable (rows y)", y.rows(),
                   "Size of kernel scales (w)", w.size());
  check_size_match(function, "Size of random variable", y.cols(),
                   "rows of covariance parameter", Sigma.rows());
  check_positive(function, "Kernel rows", Sigma.rows());

  const auto& Sigma_ref = to_ref(Sigma);
  const auto& w_ref = to_ref(w);
  const auto& y_t_ref = to_ref(y.transpose());

  check_finite(function, "Kernel", Sigma_ref);
  check_symmetric(function, "Kernel", Sigma_ref);
  check_positive_finite(function, "Kernel scales", w_ref);
  check_finite(function, "Random variable", y_t_ref);

  auto ldlt_Sigma = make_ldlt_factor(Sigma_ref);
  check_ldlt_factor(function, "LDLT_Factor of Sigma", ldlt_Sigma);

  T_lp lp(0.0);

  if (y.rows() == 0) {
    return lp;
  }

  if (include_summand<propto>::value) {
    lp += NEG_LOG_SQRT_TWO_PI * y.size();
  }

  if (include_summand<propto, T_covar>::value) {
    lp -= 0.5 * log_determinant_ldlt(ldlt_Sigma) * y.rows();
  }

  if (include_summand<propto, T_w>::value) {
    lp += (0.5 * y.cols()) * sum(log(w_ref));
  }

  if (include_summand<propto, T_y, T_w, T_covar>::value) {
    lp -= 0.5 * trace_gen_inv_quad_form_ldlt(w_ref, ldlt_Sigma, y_t_ref);
  }

  return lp;
}

template <typename T_y, typename T_covar, typename T_w>
inline return_type_t<T_y, T_covar, T_w> multi_gp_lpdf(const T_y& y,
                                                      const T_covar& Sigma,
                                                      const T_w& w) {
  return multi_gp_lpdf<false>(y, Sigma, w);
}

}  // namespace math
}  // namespace stan
#endif
