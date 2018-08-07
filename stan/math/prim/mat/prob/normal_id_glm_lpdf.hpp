#ifndef STAN_MATH_PRIM_MAT_PROB_NORMAL_ID_GLM_LPDF_HPP
#define STAN_MATH_PRIM_MAT_PROB_NORMAL_ID_GLM_LPDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/mat/meta/assign_to_matrix_or_broadcast_array.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PDF of the Generalized Linear Model (GLM)
 * with Normal distribution and id link function.
 * If containers are supplied, returns the log sum of the probabilities.
 * @tparam T_n type of vector of dependent variables (labels);
 * this can also be a single value;
 * @tparam T_x type of the matrix of independent variables (features); this
 * should be an Eigen::Matrix type whose number of rows should match the
 * length of n and whose number of columns should match the length of beta
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as n) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @tparam T_scale type of the (positive) scale(s);
 * this can be a vector (of the same length as n, for heteroskedasticity)
 * or a scalar.
 * @param n vector parameter
 * @param x design matrix
 * @param beta weight vector
 * @param alpha intercept (in log odds)
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if the scale is not positive.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_n, typename T_x, typename T_alpha,
          typename T_beta, typename T_scale>
typename return_type<T_n, T_x, T_alpha, T_beta, T_scale>::type
normal_id_glm_lpdf(const T_n &n, const T_x &x, const T_alpha &alpha,
                   const T_beta &beta, const T_scale &sigma) {
  static const char *function = "normal_id_glm_lpdf";
  typedef typename stan::partials_return_type<T_n, T_x, T_alpha, T_beta,
                                              T_scale>::type T_partials_return;

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;

  if (!(stan::length(n) && stan::length(x) && stan::length(beta)
        && stan::length(sigma)))
    return 0.0;

  T_partials_return logp(0.0);

  check_finite(function, "Vector of dependent variables", n);
  check_finite(function, "Matrix of independent variables", x);
  check_finite(function, "Weight vector", beta);
  check_finite(function, "Intercept", alpha);
  check_positive(function, "Scale vector", sigma);
  check_consistent_sizes(function, "Rows in matrix of independent variables",
                         x.col(0), "Vector of dependent variables", n);
  if (is_vector<T_scale>::value)
    check_consistent_sizes(function, "Vector of scale parameters", sigma,
                           "Vector of dependent variables", n);
  check_consistent_sizes(function, "Columns in matrix of independent variables",
                         x.row(0), "Weight vector", beta);
  if (is_vector<T_alpha>::value)
    check_consistent_sizes(function, "Vector of intercepts", alpha,
                           "Vector of dependent variables", n);

  if (!include_summand<propto, T_n, T_x, T_alpha, T_beta, T_scale>::value)
    return 0.0;

  const size_t N = x.col(0).size();
  const size_t M = x.row(0).size();

  scalar_seq_view<T_scale> sigma_vec(sigma);
  Array<T_partials_return, Dynamic, 1> inv_sigma(N, 1);
  for (size_t m = 0; m < N; ++m) {
    inv_sigma[m] = 1/value_of(sigma_vec[m]);
  }
  Matrix<T_partials_return, Dynamic, 1> beta_dbl(M, 1);
  {
    scalar_seq_view<T_beta> beta_vec(beta);
    for (size_t m = 0; m < M; ++m) {
      beta_dbl[m] = value_of(beta_vec[m]);
    }
  }

  Array<T_partials_return, Dynamic, 1> mu_dbl
      = (value_of(x) * beta_dbl).array();
  scalar_seq_view<T_alpha> alpha_vec(alpha);
  for (size_t m = 0; m < N; ++m)
    mu_dbl[m] += value_of(alpha_vec[m]);
  Array<T_partials_return, Dynamic, 1> n_minus_mu_over_sigma
      = (value_of(n).array() - mu_dbl) * inv_sigma;
  Array<T_partials_return, Dynamic, 1> n_minus_mu_over_sigma_squared
      = n_minus_mu_over_sigma.square();

  if (include_summand<propto>::value)
    logp += NEG_LOG_SQRT_TWO_PI * N;
  if (include_summand<propto, T_scale>::value) {
    for (size_t m = 0; m < N; ++m) {
      logp -= log(value_of(sigma_vec[m]));
    }
  }
  if (include_summand<propto, T_n, T_x, T_alpha, T_beta, T_scale>::value)
    logp -= 0.5 * n_minus_mu_over_sigma_squared.sum();

  // Compute the necessary derivatives.
  operands_and_partials<T_n, T_x, T_alpha, T_beta, T_scale> ops_partials(
      n, x, alpha, beta, sigma);
  if (!(is_constant_struct<T_n>::value && is_constant_struct<T_x>::value
        && is_constant_struct<T_beta>::value
        && is_constant_struct<T_alpha>::value)) {
    Matrix<T_partials_return, Dynamic, 1> mu_derivative
        = (inv_sigma * n_minus_mu_over_sigma).matrix();
    if (!is_constant_struct<T_n>::value) {
      assign_to_matrix_or_broadcast_array(ops_partials.edge1_.partials_,
                                          -mu_derivative);
    }
    if (!is_constant_struct<T_x>::value) {
      ops_partials.edge2_.partials_ = mu_derivative * beta_dbl.transpose();
    }
    if (!is_constant_struct<T_beta>::value) {
      assign_to_matrix_or_broadcast_array(
          ops_partials.edge4_.partials_,
          value_of(x).transpose() * mu_derivative);
    }
    if (!is_constant_struct<T_alpha>::value) {
      if (is_vector<T_alpha>::value)
        assign_to_matrix_or_broadcast_array(ops_partials.edge3_.partials_,
                                            mu_derivative);
      else
        ops_partials.edge3_.partials_[0] = mu_derivative.sum();
    }
    if (!is_constant_struct<T_scale>::value) {
      if (is_vector<T_scale>::value) {
        assign_to_matrix_or_broadcast_array(
            ops_partials.edge5_.partials_,
            ((n_minus_mu_over_sigma_squared
              - Array<double, Dynamic, 1>::Ones(N, 1))
             * inv_sigma)
                .matrix());
      } else {
        ops_partials.edge5_.partials_[0] =
          ((n_minus_mu_over_sigma_squared
            - Array<double, Dynamic, 1>::Ones(N, 1))
           * inv_sigma)
              .sum();
      }
    }
  }

  return ops_partials.build(logp);
}

template <typename T_n, typename T_x, typename T_alpha, typename T_beta,
          typename T_scale>
inline typename return_type<T_n, T_x, T_alpha, T_beta, T_scale>::type
normal_id_glm_lpdf(const T_n &n, const T_x &x, const T_alpha &alpha,
                   const T_beta &beta, const T_scale &sigma) {
  return normal_id_glm_lpdf<false>(n, x, alpha, beta, sigma);
}
}  // namespace math
}  // namespace stan
#endif
