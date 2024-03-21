#ifndef STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/prob/multi_normal_log.hpp>
#include <cmath>
#include <cstdlib>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the multivariate student t density for the given y, mu,
 * nu, and scale matrix.
 *
 * This version of the function is vectorized on y and mu.
 *
 * @param y scalar vector of random variates
 * @param nu scalar degrees of freedom
 * @param mu location vector
 * @param Sigma scale matrix
 * @return The log of the multivariate student t density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of scalar.
 * @tparam T_loc Type of location.
 * @tparam T_scale Type of scale.
 * @return log probability of the multivariate student t distribution.
 */
template <bool propto, typename T_y, typename T_dof, typename T_loc,
          typename T_scale>
return_type_t<T_y, T_dof, T_loc, T_scale> multi_student_t_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& Sigma) {
  using T_scale_elem = typename scalar_type<T_scale>::type;
  using lp_type = return_type_t<T_y, T_dof, T_loc, T_scale>;
  using Eigen::Matrix;
  using std::log;
  using std::vector;
  static constexpr const char* function = "multi_student_t";
  check_not_nan(function, "Degrees of freedom parameter", nu);
  check_positive(function, "Degrees of freedom parameter", nu);
  check_finite(function, "Degrees of freedom parameter", nu);

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);

  vector_seq_view<T_y> y_vec(y);
  vector_seq_view<T_loc> mu_vec(mu);
  size_t size_vec = max_size_mvt(y, mu);
  if (size_vec == 0) {
    return 0;
  }

  int num_dims = y_vec[0].size();
  if (num_dims == 0) {
    return 0;
  }

  for (size_t i = 1, size_mvt_y = size_mvt(y); i < size_mvt_y; i++) {
    check_size_match(
        function, "Size of one of the vectors of the random variable",
        y_vec[i].size(), "Size of another vector of the random variable",
        y_vec[i - 1].size());
  }

  for (size_t i = 1, size_mvt_mu = size_mvt(mu); i < size_mvt_mu; i++) {
    check_size_match(function,
                     "Size of one of the vectors "
                     "of the location variable",
                     mu_vec[i].size(),
                     "Size of another vector of "
                     "the location variable",
                     mu_vec[i - 1].size());
  }

  check_size_match(function, "Size of random variable", num_dims,
                   "size of location parameter", mu_vec[0].size());
  check_size_match(function, "Size of random variable", num_dims,
                   "rows of scale parameter", Sigma.rows());
  check_size_match(function, "Size of random variable", num_dims,
                   "columns of scale parameter", Sigma.cols());

  for (size_t i = 0; i < size_vec; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
    check_not_nan(function, "Random variable", y_vec[i]);
  }
  const auto& Sigma_ref = to_ref(Sigma);
  check_symmetric(function, "Scale parameter", Sigma_ref);

  auto ldlt_Sigma = make_ldlt_factor(Sigma_ref);
  check_ldlt_factor(function, "LDLT_Factor of scale parameter", ldlt_Sigma);

  lp_type lp(0);

  if (include_summand<propto, T_dof>::value) {
    lp += lgamma(0.5 * (nu + num_dims)) * size_vec;
    lp -= lgamma(0.5 * nu) * size_vec;
    lp -= (0.5 * num_dims) * log(nu) * size_vec;
  }

  if (include_summand<propto>::value) {
    lp -= (0.5 * num_dims) * LOG_PI * size_vec;
  }

  using Eigen::Array;

  if (include_summand<propto, T_scale_elem>::value) {
    lp -= 0.5 * log_determinant_ldlt(ldlt_Sigma) * size_vec;
  }

  if (include_summand<propto, T_y, T_dof, T_loc, T_scale_elem>::value) {
    lp_type sum_lp_vec(0.0);
    for (size_t i = 0; i < size_vec; i++) {
      const auto& y_col = as_column_vector_or_scalar(y_vec[i]);
      const auto& mu_col = as_column_vector_or_scalar(mu_vec[i]);
      sum_lp_vec
          += log1p(trace_inv_quad_form_ldlt(ldlt_Sigma, y_col - mu_col) / nu);
    }
    lp -= 0.5 * (nu + num_dims) * sum_lp_vec;
  }
  return lp;
}

template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_dof, T_loc, T_scale> multi_student_t_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& Sigma) {
  return multi_student_t_lpdf<false>(y, nu, mu, Sigma);
}

}  // namespace math
}  // namespace stan
#endif
