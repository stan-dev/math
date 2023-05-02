#ifndef STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_NEG_BINOMIAL_2_LOG_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/digamma.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/multiply_log.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Negative-Binomial-2 distribution and log link function.
 * The idea is that neg_binomial_2_log_glm_lpmf(y, x, alpha, beta, phi) should
 * compute a more efficient version of
 * neg_binomial_2_log_lpmf(y, alpha + x * beta, phi) by using analytically
 * simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_y type of positive int vector of variates (labels);
 * this can also be a single positive integer value;
 * @tparam T_x_scalar type of a scalar in the matrix of independent variables
 * (features)
 * @tparam T_x_rows compile-time number of rows of `x`. It can be either
 * `Eigen::Dynamic` or 1.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a scalar;
 * @tparam T_precision type of the (positive) precision(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 *
 * @param y failures count scalar or vector parameter. If it is a scalar it will
 * be broadcast - used for all instances.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
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
          typename T_beta, typename T_precision,
          require_matrix_t<T_x>* = nullptr>
return_type_t<T_x, T_alpha, T_beta, T_precision> neg_binomial_2_log_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta,
    const T_precision& phi) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::exp;
  using Eigen::log1p;
  using Eigen::Matrix;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_partials_return
      = partials_return_t<T_y, T_x, T_alpha, T_beta, T_precision>;
  using T_precision_val = typename std::conditional_t<
      is_vector<T_precision>::value,
      Eigen::Array<partials_return_t<T_precision>, -1, 1>,
      partials_return_t<T_precision>>;
  using T_sum_val = typename std::conditional_t<
      is_vector<T_y>::value || is_vector<T_precision>::value,
      Eigen::Array<partials_return_t<T_y, T_precision>, -1, 1>,
      partials_return_t<T_y, T_precision>>;
  using T_theta_tmp =
      typename std::conditional_t<T_x_rows == 1, T_partials_return,
                                  Array<T_partials_return, Dynamic, 1>>;
  using T_x_ref = ref_type_if_t<!is_constant<T_x>::value, T_x>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_alpha>::value, T_alpha>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;
  using T_phi_ref
      = ref_type_if_t<!is_constant<T_precision>::value, T_precision>;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();

  static const char* function = "neg_binomial_2_log_glm_lpmf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of precision parameters", phi,
                        N_instances);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  const auto& beta_val = value_of_rec(beta_ref);
  const auto& alpha_val = value_of_rec(alpha_ref);
  const auto& beta_val_vec = to_ref(as_column_vector_or_scalar(beta_val));
  const auto& alpha_val_vec = to_ref(as_column_vector_or_scalar(alpha_val));
  check_finite(function, "Weight vector", beta_val_vec);
  check_finite(function, "Intercept", alpha_val_vec);

  if (size_zero(y, phi)) {
    return 0;
  }

  const auto& y_ref = to_ref(y);
  T_phi_ref phi_ref = phi;

  const auto& y_val = value_of_rec(y_ref);
  const auto& phi_val = value_of_rec(phi_ref);

  const auto& y_val_vec = to_ref(as_column_vector_or_scalar(y_val));
  const auto& phi_val_vec = to_ref(as_column_vector_or_scalar(phi_val));
  check_nonnegative(function, "Failures variables", y_val_vec);
  check_positive_finite(function, "Precision parameter", phi_val_vec);

  if (!include_summand<propto, T_x, T_alpha, T_beta, T_precision>::value) {
    return 0;
  }

  T_x_ref x_ref = x;

  const auto& x_val
      = to_ref_if<!is_constant<T_beta>::value>(value_of_rec(x_ref));

  const auto& y_arr = as_array_or_scalar(y_val_vec);
  const auto& phi_arr = as_array_or_scalar(phi_val_vec);

  Array<T_partials_return, Dynamic, 1> theta(N_instances);
  if (T_x_rows == 1) {
    T_theta_tmp theta_tmp
        = forward_as<T_theta_tmp>((x_val * beta_val_vec)(0, 0));
    theta = theta_tmp + as_array_or_scalar(alpha_val_vec);
  } else {
    theta = (x_val * beta_val_vec).array();
    theta += as_array_or_scalar(alpha_val_vec);
  }
  check_finite(function, "Matrix of independent variables", theta);
  T_precision_val log_phi = log(phi_arr);
  Array<T_partials_return, Dynamic, 1> logsumexp_theta_logphi
      = (theta > log_phi)
            .select(theta + log1p(exp(log_phi - theta)),
                    log_phi + log1p(exp(theta - log_phi)));

  T_sum_val y_plus_phi = y_arr + phi_arr;

  // Compute the log-density.
  T_partials_return logp(0);
  if (include_summand<propto>::value) {
    if (is_vector<T_y>::value) {
      logp -= sum(lgamma(y_arr + 1.0));
    } else {
      logp -= sum(lgamma(y_arr + 1.0)) * N_instances;
    }
  }
  if (include_summand<propto, T_precision>::value) {
    if (is_vector<T_precision>::value) {
      scalar_seq_view<decltype(phi_val_vec)> phi_vec(phi_val_vec);
      for (size_t n = 0; n < N_instances; ++n) {
        logp += multiply_log(phi_vec[n], phi_vec[n]) - lgamma(phi_vec[n]);
      }
    } else {
      logp += N_instances
              * (multiply_log(forward_as<double>(phi_val),
                              forward_as<double>(phi_val))
                 - lgamma(forward_as<double>(phi_val)));
    }
  }
  logp -= sum(y_plus_phi * logsumexp_theta_logphi);

  if (include_summand<propto, T_x, T_alpha, T_beta>::value) {
    logp += sum(y_arr * theta);
  }
  if (include_summand<propto, T_precision>::value) {
    if (is_vector<T_y>::value || is_vector<T_precision>::value) {
      logp += sum(lgamma(y_plus_phi));
    } else {
      logp += sum(lgamma(y_plus_phi)) * N_instances;
    }
  }

  // Compute the necessary derivatives.
  auto ops_partials
      = make_partials_propagator(x_ref, alpha_ref, beta_ref, phi_ref);
  if (!is_constant_all<T_x, T_beta, T_alpha, T_precision>::value) {
    Array<T_partials_return, Dynamic, 1> theta_exp = theta.exp();
    if (!is_constant_all<T_x, T_beta, T_alpha>::value) {
      Matrix<T_partials_return, Dynamic, 1> theta_derivative
          = y_arr - theta_exp * y_plus_phi / (theta_exp + phi_arr);
      if (!is_constant_all<T_beta>::value) {
        if (T_x_rows == 1) {
          edge<2>(ops_partials).partials_
              = forward_as<Matrix<T_partials_return, 1, Dynamic>>(
                  theta_derivative.sum() * x_val);
        } else {
          edge<2>(ops_partials).partials_
              = x_val.transpose() * theta_derivative;
        }
      }
      if (!is_constant_all<T_x>::value) {
        if (T_x_rows == 1) {
          edge<0>(ops_partials).partials_
              = forward_as<Array<T_partials_return, Dynamic, T_x_rows>>(
                  beta_val_vec * theta_derivative.sum());
        } else {
          edge<0>(ops_partials).partials_
              = (beta_val_vec * theta_derivative.transpose()).transpose();
        }
      }
      if (!is_constant_all<T_alpha>::value) {
        if (is_vector<T_alpha>::value) {
          partials<1>(ops_partials) = std::move(theta_derivative);
        } else {
          partials<1>(ops_partials)[0] = sum(theta_derivative);
        }
      }
    }
    if (!is_constant_all<T_precision>::value) {
      if (is_vector<T_precision>::value) {
        edge<3>(ops_partials).partials_
            = 1 - y_plus_phi / (theta_exp + phi_arr) + log_phi
              - logsumexp_theta_logphi + digamma(y_plus_phi) - digamma(phi_arr);
      } else {
        partials<3>(ops_partials)[0]
            = N_instances
              + sum(-y_plus_phi / (theta_exp + phi_arr) + log_phi
                    - logsumexp_theta_logphi + digamma(y_plus_phi)
                    - digamma(phi_arr));
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta,
          typename T_precision>
inline return_type_t<T_x, T_alpha, T_beta, T_precision>
neg_binomial_2_log_glm_lpmf(const T_y& y, const T_x& x, const T_alpha& alpha,
                            const T_beta& beta, const T_precision& phi) {
  return neg_binomial_2_log_glm_lpmf<false>(y, x, alpha, beta, phi);
}
}  // namespace math
}  // namespace stan
#endif
