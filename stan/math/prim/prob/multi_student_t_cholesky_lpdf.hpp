#ifndef STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_STUDENT_T_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/is_inf.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1p.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp>
#include <cmath>
#include <cstdlib>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the multivariate student t density for the given y, mu,
 * nu, and a Cholesky factor L of the scale matrix.
 * Sigma = LL', a square, semi-positive definite matrix.
 *
 * This version of the function is vectorized on y and mu.
 *
 * @param y scalar vector of random variates
 * @param nu scalar degrees of freedom
 * @param mu location vector
 * @param L cholesky decomposition of the scale matrix
 * of the multivariate student t distribution.
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
return_type_t<T_y, T_dof, T_loc, T_scale> multi_student_t_cholesky_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& L) {
  using T_scale_elem = typename scalar_type<T_scale>::type;
  using Eigen::Matrix;
  using std::log;
  using std::vector;
  using T_partials_return = partials_return_t<T_y, T_dof, T_loc, T_covar>;
  using matrix_partials_t
      = Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using T_y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_L_ref = ref_type_t<T_scale>;

  static const char* function = "multi_student_t_cholesky";

  size_t num_y = size_mvt(y);
  size_t num_mu = size_mvt(mu);

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);

  T_y_ref y_ref = y;
  T_nu_ref nu_ref = nu;
  T_mu_ref mu_ref = mu;
  T_L_ref L_ref = L;
  vector_seq_view<T_y_ref> y_vec(y_ref);
  vector_seq_view<T_mu_ref> mu_vec(mu_ref);
  const size_t size_vec = max_size_mvt(y_ref, mu_ref);

  if (size_vec == 0) {
    return lp_type(0);
  }

  check_not_nan(function, "Degrees of freedom parameter", nu_ref);
  check_positive(function, "Degrees of freedom parameter", nu_ref);
  check_finite(function, "Degrees of freedom parameter", nu_ref);

  const int num_dims = y_vec[0].size();
  const int size_mu = mu_vec[0].size();

  for (size_t i = 1, size_mvt_y = num_y; i < size_mvt_y; i++) {
    check_size_match(
        function, "Size of one of the vectors of the random variable",
        y_vec[i].size(), "Size of the first vector of the random variable",
        num_dims);
  }

  for (size_t i = 1, size_mvt_mu = num_mu; i < size_mvt_mu; i++) {
    check_size_match(function,
                     "Size of one of the vectors "
                     "of the location variable",
                     mu_vec[i].size(),
                     "Size of another vector of "
                     "the location variable",
                     size_mu);
  }

  check_size_match(function, "Size of random variable", mu_vec[0].size(),
                   "rows of scale parameter", L.rows());
  check_size_match(function, "Size of random variable", num_dims,
                   "size of location parameter", mu_vec[0].size());
  check_size_match(function, "Size of random variable", num_dims,
                   "rows of scale parameter", L.rows());
  check_size_match(function, "Size of random variable", num_dims,
                   "columns of scale parameter", L.cols());

  if (unlikely(num_dims == 0)) {
    return lp_type(0);
  }

  for (size_t i = 0; i < size_vec; i++) {
    check_finite(function, "Location parameter", mu_vec[i]);
    check_not_nan(function, "Random variable", y_vec[i]);
  }

  check_cholesky_factor(function, "scale parameter", L_ref);

  T_partials_return lp(0);

  operands_and_partials<T_y_ref, T_dof_ref, T_mu_ref, T_L_ref> ops_partials(y_ref, nu_ref, mu_ref,
                                                                 L_ref);
  if (include_summand<propto>::value) {
    lp += -(0.5 * num_dims) * LOG_PI * size_vec;
  }

if (include_summand<propto, T_y, T_loc, T_covar_elem>::value) {

  Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>
        y_val_minus_mu_val(num_dims, size_vec);

  for (size_t i = 0; i < num_dims; i++) {
      decltype(auto) y_val = as_value_column_vector_or_scalar(y_vec[i]);
      decltype(auto) mu_val = as_value_column_vector_or_scalar(mu_vec[i]);
      y_val_minus_mu_val.col(i) = y_val - mu_val;
  }

  matrix_partials_t half;
  matrix_partials_t scaled_diff;
  auto nu_val = value_of(nu);

  // If the covariance is not autodiff, we can avoid computing a matrix
  // inverse
  if (is_constant<T_covar_elem>::value) {
    matrix_partials_t L_val = value_of(L_ref);

    half = mdivide_left_tri<Eigen::Lower>(L_val, y_val_minus_mu_val)
                 .transpose();

    scaled_diff = mdivide_right_tri<Eigen::Lower>(half, L_val).transpose();

    if (include_summand<propto>::value) {
        lp += -sum(log(L_val.diagonal())) * size_vec;
    }
  } else {
     matrix_partials_t inv_L_val
          = mdivide_left_tri<Eigen::Lower>(value_of(L_ref));
     matrix_partials_t scale_vector;

      half = (inv_L_val.template triangularView<Eigen::Lower>()
              * y_val_minus_mu_val)
                 .transpose();

      scaled_diff = (half * inv_L_val.template triangularView<Eigen::Lower>())
                        .transpose();

      scale_vector = (1 + num_dims / nu_val) / (1 + half.transpose() * half / nu_val);
      
      lp += sum(log(inv_L_val.diagonal())) * size_vec;

      ops_partials.edge3_.partials_ += -size_vec * inv_L_val.transpose();
  }

  if (include_summand<propto, T_dof>::value) {
    lp += lgamma(0.5 * (nu + num_dims)) * size_vec;
    lp += -lgamma(0.5 * nu) * size_vec;
    lp += -(0.5 * num_dims) * log(nu) * size_vec;
  }

  for (size_t i = 0; i < size_vec; i++) {
      const auto& y_col = as_column_vector_or_scalar(y_vec[i]);
      const auto& mu_col = as_column_vector_or_scalar(mu_vec[i]);
      sum_lp_vec += log1p(
          dot_self(mdivide_left_tri<Eigen::Lower>(L_ref, y_col - mu_col)) / nu);
  }
  
    lp += -0.5 * (nu + num_dims) * log1p(columns_dot_self(half) / nu);
  
  return lp;
    }
}

template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
inline return_type_t<T_y, T_dof, T_loc, T_scale> multi_student_t_cholesky_lpdf(
    const T_y& y, const T_dof& nu, const T_loc& mu, const T_scale& L) {
  return multi_student_t_cholesky_lpdf<false>(y, nu, mu, L);
}

}  // namespace math
}  // namespace stan
#endif
