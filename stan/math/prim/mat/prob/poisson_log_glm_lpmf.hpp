#ifndef STAN_MATH_PRIM_MAT_PROB_POISSON_LOG_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_POISSON_LOG_GLM_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Poisson distribution and log link function.
 * The idea is that poisson_log_glm_lpmf(y, x, alpha, beta) should
 * compute a more efficient version of poisson_log_lpmf(y, alpha + x * beta)
 * by using analytically simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 * @tparam T_y type of vector of variates (labels), integers >=0;
 * this can also be a single positive integer;
 * @tparam T_x type of the matrix of covariates (features); this
 * should be an Eigen::Matrix type whose number of rows should match the
 * length of y and whose number of columns should match the length of beta
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @param y positive integer vector parameter
 * @param x design matrix
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is negative.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta>
typename return_type<T_x, T_alpha, T_beta>::type poisson_log_glm_lpmf(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta) {
  static const char *function = "poisson_log_glm_lpmf";
  typedef typename stan::partials_return_type<T_y, T_x, T_alpha, T_beta>::type
      T_partials_return;

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;

  if (!(stan::length(y) && stan::length(x) && stan::length(beta)))
    return 0.0;

  T_partials_return logp(0.0);

  const size_t N = x.col(0).size();
  const size_t M = x.row(0).size();

  check_nonnegative(function, "Vector of dependent variables", y);
  check_finite(function, "Matrix of independent variables", x);
  check_finite(function, "Weight vector", beta);
  check_finite(function, "Intercept", alpha);
  check_consistent_size(function, "Vector of dependent variables", y, N);
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_alpha>::value)
    check_consistent_sizes(function, "Vector of intercepts", alpha,
                           "Vector of dependent variables", y);

  if (!include_summand<propto, T_x, T_alpha, T_beta>::value)
    return 0.0;

  Matrix<T_partials_return, Dynamic, 1> beta_dbl(M, 1);
  {
    scalar_seq_view<T_beta> beta_vec(beta);
    for (size_t m = 0; m < M; ++m) {
      beta_dbl[m] = value_of(beta_vec[m]);
    }
  }

  Matrix<T_partials_return, Dynamic, 1> theta_dbl = value_of(x) * beta_dbl;

  Matrix<T_partials_return, Dynamic, 1> theta_derivative(N);
  T_partials_return theta_derivative_sum = 0;
  scalar_seq_view<T_alpha> alpha_vec(alpha);
  scalar_seq_view<T_y> y_seq_view(y);
  for (size_t n = 0; n < N; ++n) {
    // Compute the derivative wrt theta.
    double theta_dbl_n = theta_dbl[n] + value_of(alpha_vec[n]);
    double exp_theta_n = exp(theta_dbl_n);
    theta_derivative[n] = y_seq_view[n] - exp_theta_n;
    if (!is_vector<T_alpha>::value)
      theta_derivative_sum += theta_derivative[n];

    // Compute the log-density.
    if (!(theta_dbl_n == -std::numeric_limits<double>::infinity()
          && y_seq_view[n] == 0)) {
      if (include_summand<propto>::value)
        logp -= lgamma(y_seq_view[n] + 1.0);
      if (include_summand<propto, T_partials_return>::value)
        logp += y_seq_view[n] * theta_dbl_n - exp_theta_n;
    }
  }

  // Compute the necessary derivatives.
  operands_and_partials<T_x, T_alpha, T_beta> ops_partials(x, alpha, beta);
  if (!(is_constant_struct<T_x>::value && is_constant_struct<T_beta>::value
        && is_constant_struct<T_alpha>::value)) {
    if (!is_constant_struct<T_beta>::value) {
      ops_partials.edge3_.partials_
          = value_of(x).transpose() * theta_derivative;
    }
    if (!is_constant_struct<T_x>::value) {
      ops_partials.edge1_.partials_ = theta_derivative * beta_dbl.transpose();
    }
    if (!is_constant_struct<T_alpha>::value) {
      if (is_vector<T_alpha>::value)
        ops_partials.edge2_.partials_ = theta_derivative;
      else
        ops_partials.edge2_.partials_[0] = theta_derivative_sum;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta>
inline typename return_type<T_x, T_alpha, T_beta>::type poisson_log_glm_lpmf(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta) {
  return poisson_log_glm_lpmf<false>(y, x, alpha, beta);
}
}  // namespace math
}  // namespace stan
#endif
