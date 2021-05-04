#ifndef STAN_MATH_PRIM_PROB_LKJ_COV_LPDF_HPP
#define STAN_MATH_PRIM_PROB_LKJ_COV_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/prob/lognormal_lpdf.hpp>
#include <stan/math/prim/prob/lkj_corr_lpdf.hpp>

namespace stan {
namespace math {

/**
  @deprecated
  */
// LKJ_cov(y|mu, sigma, eta) [ y covariance matrix (not correlation matrix)
//                         mu vector, sigma > 0 vector, eta > 0 ]
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape, require_eigen_matrix_dynamic_t<T_y>* = nullptr,
          require_all_eigen_col_vector_t<T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_shape> lkj_cov_lpdf(const T_y& y,
                                                         const T_loc& mu,
                                                         const T_scale& sigma,
                                                         const T_shape& eta) {
  static const char* function = "lkj_cov_lpdf";
  check_size_match(function, "Rows of location parameter", mu.rows(),
                   "columns of scale parameter", sigma.rows());
  check_square(function, "random variable", y);
  check_size_match(function, "Rows of random variable", y.rows(),
                   "rows of location parameter", mu.rows());
  const auto& y_ref = to_ref(y);
  const auto& mu_ref = to_ref(mu);
  const auto& sigma_ref = to_ref(sigma);
  check_positive(function, "Shape parameter", eta);
  check_finite(function, "Location parameter", mu_ref);
  check_finite(function, "Scale parameter", sigma_ref);
  check_finite(function, "Covariance matrix", y_ref);

  return_type_t<T_y, T_loc, T_scale, T_shape> lp(0.0);

  const unsigned int K = y.rows();
  const Eigen::Array<value_type_t<T_y>, Eigen::Dynamic, 1> sds
      = y_ref.diagonal().array().sqrt();
  for (unsigned int k = 0; k < K; k++) {
    lp += lognormal_lpdf<propto>(sds(k), mu_ref(k), sigma_ref(k));
  }
  if (stan::is_constant_all<T_shape>::value && eta == 1.0) {
    // no need to rescale y into a correlation matrix
    lp += lkj_corr_lpdf<propto>(y_ref, eta);
    return lp;
  }
  Eigen::DiagonalMatrix<value_type_t<T_y>, Eigen::Dynamic> D(K);
  D.diagonal() = sds.inverse();
  lp += lkj_corr_lpdf<propto>(D * y_ref * D, eta);
  return lp;
}

/**
  @deprecated
  */
// LKJ_Cov(y|mu, sigma, eta) [ y covariance matrix (not correlation matrix)
//                         mu scalar, sigma > 0 scalar, eta > 0 ]
template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_shape, require_eigen_matrix_dynamic_t<T_y>* = nullptr,
          require_all_stan_scalar_t<T_loc, T_scale>* = nullptr>
return_type_t<T_y, T_loc, T_scale, T_shape> lkj_cov_lpdf(const T_y& y,
                                                         const T_loc& mu,
                                                         const T_scale& sigma,
                                                         const T_shape& eta) {
  static const char* function = "lkj_cov_lpdf";
  check_positive(function, "Shape parameter", eta);
  check_finite(function, "Location parameter", mu);
  check_finite(function, "Scale parameter", sigma);
  const auto& y_ref = to_ref(y);
  check_finite(function, "Covariance matrix", y_ref);

  return_type_t<T_y, T_loc, T_scale, T_shape> lp(0.0);

  const unsigned int K = y.rows();
  const Eigen::Array<value_type_t<T_y>, Eigen::Dynamic, 1> sds
      = y_ref.diagonal().array().sqrt();
  for (unsigned int k = 0; k < K; k++) {
    lp += lognormal_lpdf<propto>(sds(k), mu, sigma);
  }
  if (stan::is_constant_all<T_shape>::value && eta == 1.0) {
    // no need to rescale y into a correlation matrix
    lp += lkj_corr_lpdf<propto>(y_ref, eta);
    return lp;
  }
  Eigen::DiagonalMatrix<value_type_t<T_y>, Eigen::Dynamic> D(K);
  D.diagonal() = sds.inverse();
  lp += lkj_corr_lpdf<propto>(D * y_ref * D, eta);
  return lp;
}

/**
  @deprecated
  */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> lkj_cov_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& eta) {
  return lkj_cov_lpdf<false>(y, mu, sigma, eta);
}

}  // namespace math
}  // namespace stan
#endif
