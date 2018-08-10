#ifndef STAN_MATH_PRIM_MAT_PROB_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/mat/meta/is_vector.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Negative-Binomial-2 distribution and log link function.
 * The idea is that neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi) should
 * compute a more efficient version of
 * neg_binomial_2_log_lpmf(y, alpha + x * beta, phi) by using analytically
 * simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 * @tparam T_y type of positive int vector of variates (labels);
 * this can also be a single positive integer value;
 * @tparam T_x type of the matrix of covariates (features); this
 * should be an Eigen::Matrix type whose number of rows should match the
 * length of y and whose number of columns should match the length of beta
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a scalar;
 * @tparam T_precision type of the (positive) precision(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 * @param y failures count vector parameter
 * @param x design matrix
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param phi (vector of) precision parameter(s)
 * @return log probability or log sum of probabilities
 * @throw std::invalid_argument if container sizes mismatch.
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if phi is infinite or non-positive.
 * @throw std::domain_error if y is negative.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, typename T_precision>
typename return_type<T_x, T_alpha, T_beta, T_precision>::type
neg_binomial_2_log_glm_lpmf(const T_y& y, const T_x& x, const T_alpha& alpha,
                            const T_beta& beta, const T_precision& phi) {
  static const char* function = "neg_binomial_2_log_glm_lpmf";
  typedef
      typename stan::partials_return_type<T_y, T_x, T_alpha, T_beta,
                                          T_precision>::type T_partials_return;

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;

  if (!(stan::length(y) && stan::length(x) && stan::length(beta)
        && stan::length(phi)))
    return 0.0;

  T_partials_return logp(0.0);

  const size_t N = x.col(0).size();
  const size_t M = x.row(0).size();

  check_nonnegative(function, "Failures variables", y);
  check_finite(function, "Matrix of independent variables", x);
  check_finite(function, "Weight vector", beta);
  check_finite(function, "Intercept", alpha);
  check_positive_finite(function, "Precision parameter", phi);
  check_consistent_size(function, "Vector of dependent variables", y, N);
  check_consistent_size(function, "Weight vector", beta, M);
  if (is_vector<T_precision>::value)
    check_consistent_sizes(function, "Vector of precision parameters", phi,
                           "Vector of dependent variables", y);
  if (is_vector<T_alpha>::value)
    check_consistent_sizes(function, "Vector of intercepts", alpha,
                           "Vector of dependent variables", y);

  if (!include_summand<propto, T_x, T_alpha, T_beta, T_precision>::value)
    return 0.0;

  Array<T_partials_return, Dynamic, 1> y_arr(N, 1);
  Array<T_partials_return, Dynamic, 1> phi_arr(N, 1);
  {
    scalar_seq_view<T_y> y_vec(y);
    scalar_seq_view<T_precision> phi_vec(phi);
    for (size_t n = 0; n < N; ++n) {
      y_arr[n] = y_vec[n];
      phi_arr[n] = value_of(phi_vec[n]);
    }
  }

  Matrix<T_partials_return, Dynamic, 1> beta_dbl(M, 1);
  {
    scalar_seq_view<T_beta> beta_vec(beta);
    for (size_t m = 0; m < M; ++m) {
      beta_dbl[m] = value_of(beta_vec[m]);
    }
  }

  Array<T_partials_return, Dynamic, 1> theta_dbl
      = (value_of(x) * beta_dbl).array();
  scalar_seq_view<T_alpha> alpha_vec(alpha);
  for (size_t n = 0; n < N; ++n)
    theta_dbl[n] += value_of(alpha_vec[n]);
  Array<T_partials_return, Dynamic, 1> log_phi = phi_arr.log();
  Array<T_partials_return, Dynamic, 1> logsumexp_eta_logphi
      = theta_dbl.binaryExpr(log_phi, [](const T_partials_return& xx,
                                         const T_partials_return& yy) {
          return log_sum_exp(xx, yy);
        });
  Array<T_partials_return, Dynamic, 1> y_plus_phi = y_arr + phi_arr;

  // Compute the log-density.
  if (include_summand<propto>::value) {
    logp -= (y_arr + Array<double, Dynamic, 1>::Ones(N, 1))
                .unaryExpr(
                    [](const T_partials_return& xx) { return lgamma(xx); })
                .sum();
  }
  if (include_summand<propto, T_precision>::value) {
    for (size_t n = 0; n < N; ++n)
      logp += multiply_log(phi_arr[n], phi_arr[n]) - lgamma(phi_arr[n]);
  }
  if (include_summand<propto, T_x, T_alpha, T_beta, T_precision>::value)
    logp -= (y_plus_phi * logsumexp_eta_logphi).sum();
  if (include_summand<propto, T_x, T_alpha, T_beta>::value)
    logp += (y_arr * theta_dbl).sum();
  if (include_summand<propto, T_precision>::value) {
    logp += y_plus_phi
                .unaryExpr(
                    [](const T_partials_return& xx) { return lgamma(xx); })
                .sum();
  }

  // Compute the necessary derivatives.
  operands_and_partials<T_x, T_alpha, T_beta, T_precision> ops_partials(
      x, alpha, beta, phi);

  if (!(is_constant_struct<T_x>::value && is_constant_struct<T_beta>::value
        && is_constant_struct<T_alpha>::value)) {
    Matrix<T_partials_return, Dynamic, 1> theta_derivative(N, 1);
    theta_derivative = (y_arr
                        - (y_plus_phi
                           / (phi_arr / (theta_dbl.exp())
                              + Array<double, Dynamic, 1>::Ones(N, 1))))
                           .matrix();
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
        ops_partials.edge2_.partials_[0] = theta_derivative.sum();
    }
  }
  if (!is_constant_struct<T_precision>::value) {
    if (is_vector<T_precision>::value) {
      ops_partials.edge4_.partials_
          = (Array<double, Dynamic, 1>::Ones(N, 1)
             - y_plus_phi / (theta_dbl.exp() + phi_arr) + log_phi
             - logsumexp_eta_logphi
             - phi_arr.unaryExpr(
                   [](const T_partials_return& xx) { return digamma(xx); })
             + y_plus_phi.unaryExpr(
                   [](const T_partials_return& xx) { return digamma(xx); }))
                .matrix();
    } else {
      ops_partials.edge4_.partials_[0]
          = (Array<double, Dynamic, 1>::Ones(N, 1)
             - y_plus_phi / (theta_dbl.exp() + phi_arr) + log_phi
             - logsumexp_eta_logphi
             - phi_arr.unaryExpr(
                   [](const T_partials_return& xx) { return digamma(xx); })
             + y_plus_phi.unaryExpr(
                   [](const T_partials_return& xx) { return digamma(xx); }))
                .sum();
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta,
          typename T_precision>
inline typename return_type<T_x, T_alpha, T_beta, T_precision>::type
neg_binomial_2_log_glm_lpmf(const T_y& y, const T_x& x, const T_alpha& alpha,
                            const T_beta& beta, const T_precision& phi) {
  return neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta, phi);
}
}  // namespace math
}  // namespace stan
#endif
