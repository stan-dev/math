#ifndef STAN_MATH_PRIM_PROB_BINOMIAL_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BINOMIAL_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/max.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_vector_or_scalar.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Binomial distribution and logit link function.
 * The idea is that binomial_logit_glm_lpmf(n | N, x, alpha, beta) should
 * compute a more efficient version of
 * binomial_logit_lpmf(y | N, alpha + x * beta) by using analytically
 * simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_n type of binary vector of successes variables;
 * this can also be a single binary value;
 * @tparam T_N type of binary vector of population size variables;
 * this can also be a single binary value;
 * @tparam T_x type of the matrix of independent variables (features)
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector
 *
 * @param n binary scalar or vector parameter. If it is a scalar it will be
 * broadcast - used for all instances.
 * @param N binary scalar or vector parameter. If it is a scalar it will be
 * broadcast - used for all instances.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if n is negative or greater than N
 * @throw std::domain_error if N is negative
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_n, typename T_N, typename T_x,
          typename T_alpha, typename T_beta, require_matrix_t<T_x>* = nullptr>
return_type_t<T_x, T_alpha, T_beta> binomial_logit_glm_lpmf(
    const T_n& n, const T_N& N, const T_x& x, const T_alpha& alpha,
    const T_beta& beta) {
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_xbeta_partials = partials_return_t<T_x, T_beta>;
  using T_partials_return = partials_return_t<T_x, T_alpha, T_beta>;
  using T_xbeta_tmp =
      typename std::conditional_t<T_x_rows == 1, T_xbeta_partials,
                                  Eigen::Array<T_xbeta_partials, -1, 1>>;
  using T_n_ref = ref_type_if_t<!is_constant<T_n>::value, T_n>;
  using T_N_ref = ref_type_if_t<!is_constant<T_N>::value, T_N>;
  using T_x_ref = ref_type_if_t<!is_constant<T_x>::value, T_x>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_alpha>::value, T_alpha>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;

  T_n_ref n_ref = n;
  T_N_ref N_ref = N;
  T_x_ref x_ref = x;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  if (size_zero(n, N, alpha, beta, x)) {
    return 0;
  }

  if (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
    return 0;
  }

  const size_t N_instances = max_size(n, N, x.col(0), alpha);
  const size_t N_attributes = x.cols();

  static const char* function = "binomial_logit_glm_lpmf";
  check_consistent_sizes(function, "Successes variable", n,
                         "Population size parameter", N);
  check_consistent_size(function, "Successes variable", n, N_instances);
  check_consistent_size(function, "Population size parameter", N, N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);

  auto&& n_val = as_value_column_array_or_scalar(n_ref);
  auto&& N_val = as_value_column_array_or_scalar(N_ref);

  check_bounded(function, "Successes variable", n_val, 0, N_val);
  check_nonnegative(function, "Population size parameter", N_val);

  auto&& alpha_val = as_value_column_array_or_scalar(alpha_ref);
  auto&& beta_val = as_value_column_vector_or_scalar(beta_ref);
  auto&& x_val = value_of(x_ref);
  Eigen::Array<T_partials_return, -1, 1> theta(N_instances);
  if (T_x_rows == 1) {
    theta = forward_as<T_xbeta_tmp>((x_val * beta_val)(0, 0)) + alpha_val;
  } else {
    theta = (x_val * beta_val).array() + alpha_val;
  }

  constexpr bool gradients_calc = !is_constant_all<T_beta, T_x, T_alpha>::value;
  auto&& log_inv_logit_theta = to_ref_if<gradients_calc>(log_inv_logit(theta));

  T_partials_return logp = sum(n_val * log_inv_logit_theta
                               + (N_val - n_val) * log1m_inv_logit(theta));

  using std::isfinite;
  if (!isfinite(logp)) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables", x);
  }

  if (include_summand<propto, T_n, T_N>::value) {
    size_t broadcast_n = max_size(N, n) == N_instances ? 1 : N_instances;
    logp += sum(binomial_coefficient_log(N_val, n_val)) * broadcast_n;
  }

  auto ops_partials = make_partials_propagator(x_ref, alpha_ref, beta_ref);
  if (gradients_calc) {
    Eigen::Matrix<T_partials_return, -1, 1> theta_derivative
        = n_val - N_val * exp(log_inv_logit_theta);

    if (!is_constant_all<T_beta>::value) {
      if (T_x_rows == 1) {
        edge<2>(ops_partials).partials_
            = forward_as<Eigen::Matrix<T_partials_return, 1, -1>>(
                theta_derivative.sum() * x_val);
      } else {
        partials<2>(ops_partials) = x_val.transpose() * theta_derivative;
      }
    }

    if (!is_constant_all<T_x>::value) {
      if (T_x_rows == 1) {
        edge<0>(ops_partials).partials_
            = forward_as<Eigen::Array<T_partials_return, -1, T_x_rows>>(
                beta_val * theta_derivative.sum());
      } else {
        edge<0>(ops_partials).partials_
            = (beta_val * theta_derivative.transpose()).transpose();
      }
    }
    if (!is_constant_all<T_alpha>::value) {
      partials<1>(ops_partials) = theta_derivative;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_n, typename T_N, typename T_x, typename T_alpha,
          typename T_beta>
inline return_type_t<T_x, T_beta, T_alpha> binomial_logit_glm_lpmf(
    const T_n& n, const T_N& N, const T_x& x, const T_alpha& alpha,
    const T_beta& beta) {
  return binomial_logit_glm_lpmf<false>(n, N, x, alpha, beta);
}
}  // namespace math
}  // namespace stan
#endif
