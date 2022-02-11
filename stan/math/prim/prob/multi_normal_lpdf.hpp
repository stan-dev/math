#ifndef STAN_MATH_PRIM_PROB_MULTI_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_covar>
return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_covar& Sigma) {
  using T_covar_elem = typename scalar_type<T_covar>::type;
  using lp_type = return_type_t<T_y, T_loc, T_covar>;
  using Eigen::Dynamic;
  static const char* function = "multi_normal_lpdf";
  check_positive(function, "Covariance matrix rows", Sigma.rows());

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);
  size_t number_of_y = size_mvt(y);
  size_t number_of_mu = size_mvt(mu);
  if (number_of_y == 0 || number_of_mu == 0) {
    return 0.0;
  }

  lp_type lp(0.0);
  vector_seq_view<T_y> y_vec(y);
  vector_seq_view<T_loc> mu_vec(mu);
  size_t size_vec = max_size_mvt(y, mu);

  int size_y = y_vec[0].size();
  int size_mu = mu_vec[0].size();
  if (size_vec > 1) {
    for (size_t i = 1, size_mvt_y = size_mvt(y); i < size_mvt_y; i++) {
      check_size_match(function,
                       "Size of one of the vectors of "
                       "the random variable",
                       y_vec[i].size(),
                       "Size of the first vector of the "
                       "random variable",
                       size_y);
    }
    for (size_t i = 1, size_mvt_mu = size_mvt(mu); i < size_mvt_mu; i++) {
      check_size_match(function,
                       "Size of one of the vectors of "
                       "the location variable",
                       mu_vec[i].size(),
                       "Size of the first vector of the "
                       "location variable",
                       size_mu);
    }
  }

  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of covariance parameter", Sigma.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of covariance parameter", Sigma.cols());

  for (size_t i = 0; i < size_vec; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
    check_not_nan(function, "Random variable", y_vec[i]);
  }
  const auto& Sigma_ref = to_ref(Sigma);
  check_symmetric(function, "Covariance matrix", Sigma_ref);

  auto ldlt_Sigma = make_ldlt_factor(Sigma_ref);
  check_ldlt_factor(function, "LDLT_Factor of covariance parameter",
                    ldlt_Sigma);

  if (size_y == 0) {
    return lp;
  }

  if (include_summand<propto>::value) {
    lp += NEG_LOG_SQRT_TWO_PI * size_y * size_vec;
  }

  if (include_summand<propto, T_covar_elem>::value) {
    lp -= 0.5 * log_determinant_ldlt(ldlt_Sigma) * size_vec;
  }

  if (include_summand<propto, T_y, T_loc, T_covar_elem>::value) {
    lp_type sum_lp_vec(0.0);
    for (size_t i = 0; i < size_vec; i++) {
      const auto& y_col = as_column_vector_or_scalar(y_vec[i]);
      const auto& mu_col = as_column_vector_or_scalar(mu_vec[i]);
      sum_lp_vec += trace_inv_quad_form_ldlt(ldlt_Sigma, y_col - mu_col);
    }
    lp -= 0.5 * sum_lp_vec;
  }
  return lp;
}

template <typename T_y, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_covar& Sigma) {
  return multi_normal_lpdf<false>(y, mu, Sigma);
}

}  // namespace math
}  // namespace stan
#endif
