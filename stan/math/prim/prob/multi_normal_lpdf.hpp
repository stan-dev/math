#ifndef STAN_MATH_PRIM_PROB_MULTI_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MULTI_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/transpose.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_covar,
          require_any_not_vector_vt<is_stan_scalar, T_y, T_loc>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_covar>* = nullptr>
return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_covar& Sigma) {
  using T_covar_elem = typename scalar_type<T_covar>::type;
  using T_return = return_type_t<T_y, T_loc, T_covar>;
  using T_partials_return = partials_return_t<T_y, T_loc, T_covar>;
  using matrix_partials_t
      = Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using vector_partials_t = Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_Sigma_ref = ref_type_t<T_covar>;

  static const char* function = "multi_normal_lpdf";
  check_positive(function, "Covariance matrix rows", Sigma.rows());

  check_consistent_sizes_mvt(function, "y", y, "mu", mu);
  size_t number_of_y = size_mvt(y);
  size_t number_of_mu = size_mvt(mu);
  if (number_of_y == 0 || number_of_mu == 0) {
    return 0.0;
  }

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_Sigma_ref Sigma_ref = Sigma;
  vector_seq_view<T_y_ref> y_vec(y_ref);
  vector_seq_view<T_mu_ref> mu_vec(mu_ref);
  const size_t size_vec = max_size_mvt(y, mu);
  const int K = Sigma.rows();

  int size_y = y_vec[0].size();
  int size_mu = mu_vec[0].size();

  // check size consistency of all random variables y
  for (size_t i = 1, size_mvt_y = size_mvt(y); i < size_mvt_y; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the random variable",
                     y_vec[i].size(),
                     "Size of the first vector of the "
                     "random variable",
                     size_y);
  }
  // check size consistency of all means mu
  for (size_t i = 1, size_mvt_mu = size_mvt(mu); i < size_mvt_mu; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the location variable",
                     mu_vec[i].size(),
                     "Size of the first vector of the "
                     "location variable",
                     size_mu);
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
  check_symmetric(function, "Covariance matrix", Sigma_ref);

  auto ldlt_Sigma = make_ldlt_factor(value_of(Sigma_ref));

  check_ldlt_factor(function, "LDLT_Factor of covariance parameter",
                    ldlt_Sigma);

  if (unlikely(size_y == 0)) {
    return T_return(0);
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, Sigma_ref);

  T_partials_return logp(0);

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * size_y * size_vec;
  }

  if (include_summand<propto, T_y, T_loc, T_covar_elem>::value) {
    vector_partials_t half(size_vec);
    vector_partials_t y_val_minus_mu_val(size_vec);

    T_partials_return sum_lp_vec(0.0);
    for (size_t i = 0; i < size_vec; i++) {
      const auto& y_val = as_value_column_vector_or_scalar(y_vec[i]);
      const auto& mu_val = as_value_column_vector_or_scalar(mu_vec[i]);
      y_val_minus_mu_val = eval(y_val - mu_val);
      half = mdivide_left_ldlt(ldlt_Sigma, y_val_minus_mu_val);

      sum_lp_vec += dot_product(y_val_minus_mu_val, half);

      if (!is_constant_all<T_y>::value) {
        partials_vec<0>(ops_partials)[i] += -half;
      }
      if (!is_constant_all<T_loc>::value) {
        partials_vec<1>(ops_partials)[i] += half;
      }
      if (!is_constant<T_covar_elem>::value) {
        partials_vec<2>(ops_partials)[i] += 0.5 * half * half.transpose();
      }
    }

    logp += -0.5 * sum_lp_vec;

    // If the covariance is not autodiff, we can avoid computing a matrix
    // inverse
    if (is_constant<T_covar_elem>::value) {
      if (include_summand<propto>::value) {
        logp += -0.5 * log_determinant_ldlt(ldlt_Sigma) * size_vec;
      }
    } else {
      matrix_partials_t inv_Sigma
          = mdivide_left_ldlt(ldlt_Sigma, Eigen::MatrixXd::Identity(K, K));

      logp += -0.5 * log_determinant_ldlt(ldlt_Sigma) * size_vec;

      partials<2>(ops_partials) += -0.5 * size_vec * inv_Sigma;
    }
  }

  return ops_partials.build(logp);
}

template <bool propto, typename T_y, typename T_loc, typename T_covar,
          require_all_vector_vt<is_stan_scalar, T_y, T_loc>* = nullptr,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_covar>* = nullptr>
return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(const T_y& y,
                                                     const T_loc& mu,
                                                     const T_covar& Sigma) {
  using T_covar_elem = typename scalar_type<T_covar>::type;
  using T_return = return_type_t<T_y, T_loc, T_covar>;
  using T_partials_return = partials_return_t<T_y, T_loc, T_covar>;
  using matrix_partials_t
      = Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using vector_partials_t = Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using T_y_ref = ref_type_t<T_y>;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_Sigma_ref = ref_type_t<T_covar>;

  static const char* function = "multi_normal_lpdf";
  check_positive(function, "Covariance matrix rows", Sigma.rows());

  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_Sigma_ref Sigma_ref = Sigma;

  decltype(auto) y_val = as_value_column_vector_or_scalar(y_ref);
  decltype(auto) mu_val = as_value_column_vector_or_scalar(mu_ref);

  const int size_y = y_ref.size();
  const int size_mu = mu_ref.size();
  const unsigned int K = Sigma.rows();

  check_finite(function, "Location parameter", mu_val);
  check_not_nan(function, "Random variable", y_val);

  check_size_match(function, "Size of random variable", size_y,
                   "size of location parameter", size_mu);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of covariance parameter", Sigma.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of covariance parameter", Sigma.cols());

  check_symmetric(function, "Covariance matrix", Sigma_ref);

  auto ldlt_Sigma = make_ldlt_factor(value_of(Sigma_ref));
  check_ldlt_factor(function, "LDLT_Factor of covariance parameter",
                    ldlt_Sigma);

  if (unlikely(size_y == 0)) {
    return T_return(0);
  }

  auto ops_partials = make_partials_propagator(y_ref, mu_ref, Sigma_ref);

  T_partials_return logp(0);

  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * size_y;
  }

  if (include_summand<propto, T_y, T_loc, T_covar_elem>::value) {
    vector_partials_t half(size_y);
    vector_partials_t y_val_minus_mu_val = eval(y_val - mu_val);

    // If the covariance is not autodiff, we can avoid computing a matrix
    // inverse
    if (is_constant<T_covar_elem>::value) {
      half = mdivide_left_ldlt(ldlt_Sigma, y_val_minus_mu_val);

      if (include_summand<propto>::value) {
        logp += -0.5 * log_determinant_ldlt(ldlt_Sigma);
      }
    } else {
      matrix_partials_t inv_Sigma
          = mdivide_left_ldlt(ldlt_Sigma, Eigen::MatrixXd::Identity(K, K));

      half.noalias() = inv_Sigma * y_val_minus_mu_val;

      logp += -0.5 * log_determinant_ldlt(ldlt_Sigma);

      edge<2>(ops_partials).partials_
          += 0.5 * (half * half.transpose() - inv_Sigma);
    }

    logp += -0.5 * dot_product(y_val_minus_mu_val, half);

    if (!is_constant_all<T_y>::value) {
      partials<0>(ops_partials) += -half;
    }
    if (!is_constant_all<T_loc>::value) {
      partials<1>(ops_partials) += half;
    }
  }

  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_covar>
inline return_type_t<T_y, T_loc, T_covar> multi_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_covar& Sigma) {
  return multi_normal_lpdf<false>(y, mu, Sigma);
}

}  // namespace math
}  // namespace stan
#endif
