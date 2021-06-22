#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Poisson distribution and log link function.
 * The idea is that poisson_log_glm_lpmf(y, x, alpha, beta) should
 * compute a more efficient version of poisson_log_lpmf(y, alpha + x * beta)
 * by using analytically simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 *
 * @tparam T_y type of vector of variates (labels), integers >=0;
 * this can also be a single positive integer;
 * @tparam T_x type the matrix of independent variables (features)
 * @tparam T_x_rows compile-time number of rows of `x`. It can be either
 * `Eigen::Dynamic` or 1.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @param y positive integer scalar or vector parameter. If it is a scalar it
 * will be broadcast - used for all instances.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is negative.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, require_eigen_t<T_x>* = nullptr>
return_type_t<T_x, T_alpha, T_beta> poisson_log_glm_lpmf(const T_y& y,
                                                         const T_x& x,
                                                         const T_alpha& alpha,
                                                         const T_beta& beta) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::exp;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_partials_return = partials_return_t<T_y, T_x, T_alpha, T_beta>;
  using T_theta_tmp =
      typename std::conditional_t<T_x_rows == 1, T_partials_return,
                                  Array<T_partials_return, Dynamic, 1>>;
  using T_x_ref = ref_type_if_t<!is_constant<T_x>::value, T_x>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_alpha>::value, T_alpha>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();

  static const char* function = "poisson_log_glm_lpmf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);
  const auto& y_ref = to_ref(y);
  check_nonnegative(function, "Vector of dependent variables", y_ref);

  if (size_zero(y)) {
    return 0;
  }
  if (!include_summand<propto, T_x, T_alpha, T_beta>::value) {
    return 0;
  }

  T_x_ref x_ref = x;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  const auto& y_val = value_of_rec(y_ref);
  const auto& x_val
      = to_ref_if<!is_constant<T_beta>::value>(value_of_rec(x_ref));
  const auto& alpha_val = value_of_rec(alpha_ref);
  const auto& beta_val = value_of_rec(beta_ref);

  const auto& y_val_vec = to_ref(as_column_vector_or_scalar(y_val));
  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val);
  const auto& beta_val_vec = to_ref_if<!is_constant<T_x>::value>(
      as_column_vector_or_scalar(beta_val));

  Array<T_partials_return, Dynamic, 1> theta(N_instances);
  if (T_x_rows == 1) {
    T_theta_tmp theta_tmp
        = forward_as<T_theta_tmp>((x_val * beta_val_vec).coeff(0, 0));
    theta = theta_tmp + as_array_or_scalar(alpha_val_vec);
  } else {
    theta = x_val * beta_val_vec;
    theta += as_array_or_scalar(alpha_val_vec);
  }

  Matrix<T_partials_return, Dynamic, 1> theta_derivative
      = as_array_or_scalar(y_val_vec) - exp(theta.array());
  double theta_derivative_sum = sum(theta_derivative);
  if (!std::isfinite(theta_derivative_sum)) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables", theta);
  }

  T_partials_return logp(0);
  if (include_summand<propto>::value) {
    if (is_vector<T_y>::value) {
      logp -= sum(lgamma(as_array_or_scalar(y_val_vec) + 1));
    } else {
      logp -= lgamma(forward_as<double>(y_val) + 1);
    }
  }

  logp += sum(as_array_or_scalar(y_val_vec) * theta.array()
              - exp(theta.array()));

  operands_and_partials<T_x_ref, T_alpha_ref, T_beta_ref> ops_partials(
      x_ref, alpha_ref, beta_ref);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_beta>::value) {
    if (T_x_rows == 1) {
      ops_partials.edge3_.partials_
          = forward_as<Matrix<T_partials_return, 1, Dynamic>>(
              theta_derivative.sum() * x_val);
    } else {
      ops_partials.edge3_.partials_ = x_val.transpose() * theta_derivative;
    }
  }
  if (!is_constant_all<T_x>::value) {
    if (T_x_rows == 1) {
      ops_partials.edge1_.partials_
          = forward_as<Array<T_partials_return, Dynamic, T_x_rows>>(
              beta_val_vec * theta_derivative.sum());
    } else {
      ops_partials.edge1_.partials_
          = (beta_val_vec * theta_derivative.transpose()).transpose();
    }
  }
  if (!is_constant_all<T_alpha>::value) {
    if (is_vector<T_alpha>::value) {
      ops_partials.edge2_.partials_ = theta_derivative;
    } else {
      ops_partials.edge2_.partials_[0] = theta_derivative_sum;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta>
inline return_type_t<T_x, T_alpha, T_beta> poisson_log_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta) {
  return poisson_log_glm_lpmf<false>(y, x, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
