#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>

#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PMF of the Generalized Linear Model (GLM)
 * with Bernoulli distribution and logit link function.
 * The idea is that bernoulli_logit_glm_lpmf(y, x, alpha, beta) should
 * compute a more efficient version of bernoulli_logit_lpmf(y, alpha + x * beta)
 * by using analytically simplified gradients.
 * If containers are supplied, returns the log sum of the probabilities.
 * @tparam T_y type of binary vector of dependent variables (labels);
 * this can also be a single binary value;
 * @tparam T_x_scalar type of a scalar in the matrix of independent variables
 * (features)
 * @tparam T_x_rows compile-time number of rows of `x`. It can be either
 * `Eigen::Dynamic` or 1.
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @param y binary scalar or vector parameter. If it is a scalar it will be
 * broadcast - used for all instances.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if y is not binary.
 * @throw std::invalid_argument if container sizes mismatch.
 */

template <bool propto, typename T_y, typename T_x_scalar, int T_x_rows,
          typename T_alpha, typename T_beta>
return_type_t<T_x_scalar, T_alpha, T_beta> bernoulli_logit_glm_lpmf(
    const T_y &y, const Eigen::Matrix<T_x_scalar, T_x_rows, Eigen::Dynamic> &x,
    const T_alpha &alpha, const T_beta &beta) {
  static const char *function = "bernoulli_logit_glm_lpmf";

  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::log1p;
  using std::exp;

  using T_partials_return = partials_return_t<T_y, T_x_scalar, T_alpha, T_beta>;
  using T_y_val =
      typename std::conditional_t<is_vector<T_y>::value,
                                  Eigen::Matrix<partials_return_t<T_y>, -1, 1>,
                                  partials_return_t<T_y>>;
  using T_ytheta_tmp =
      typename std::conditional_t<T_x_rows == 1, T_partials_return,
                                  Array<T_partials_return, Dynamic, 1>>;

  const size_t N_instances = T_x_rows == 1 ? size(y) : x.rows();
  const size_t N_attributes = x.cols();

  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);
  check_bounded(function, "Vector of dependent variables", y, 0, 1);

  if (size_zero(y)) {
    return 0;
  }

  if (!include_summand<propto, T_x_scalar, T_alpha, T_beta>::value) {
    return 0;
  }

  T_partials_return logp(0);
  const auto &x_val = value_of_rec(x);
  const auto &y_val = value_of_rec(y);
  const auto &beta_val = value_of_rec(beta);
  const auto &alpha_val = value_of_rec(alpha);

  const auto &y_val_vec = as_column_vector_or_scalar(y_val);
  const auto &beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto &alpha_val_vec = as_column_vector_or_scalar(alpha_val);

  T_y_val signs = 2 * as_array_or_scalar(y_val_vec) - 1;

  Array<T_partials_return, Dynamic, 1> ytheta(N_instances);
  if (T_x_rows == 1) {
    T_ytheta_tmp ytheta_tmp = x_val * beta_val_vec;
    ytheta = as_array_or_scalar(signs)
             * (ytheta_tmp + as_array_or_scalar(alpha_val_vec));
  } else {
    ytheta = x_val * beta_val_vec;
    ytheta = as_array_or_scalar(signs)
             * (ytheta + as_array_or_scalar(alpha_val_vec));
  }

  // Compute the log-density and handle extreme values gracefully
  // using Taylor approximations.
  // And compute the derivatives wrt theta.
  static const double cutoff = 20.0;
  Eigen::Array<T_partials_return, Dynamic, 1> exp_m_ytheta = exp(-ytheta);
  logp += sum(
      (ytheta > cutoff)
          .select(-exp_m_ytheta,
                  (ytheta < -cutoff).select(ytheta, -log1p(exp_m_ytheta))));

  if (!std::isfinite(logp)) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables", ytheta);
  }

  operands_and_partials<Eigen::Matrix<T_x_scalar, T_x_rows, Eigen::Dynamic>,
                        T_alpha, T_beta>
      ops_partials(x, alpha, beta);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_beta, T_x_scalar, T_alpha>::value) {
    Matrix<T_partials_return, Dynamic, 1> theta_derivative
        = (ytheta > cutoff)
              .select(-exp_m_ytheta,
                      (ytheta < -cutoff)
                          .select(as_array_or_scalar(signs),
                                  as_array_or_scalar(signs) * exp_m_ytheta
                                      / (exp_m_ytheta + 1)));
    if (!is_constant_all<T_beta>::value) {
      if (T_x_rows == 1) {
        ops_partials.edge3_.partials_
            = forward_as<Matrix<T_partials_return, 1, Dynamic>>(
                theta_derivative.sum() * x_val);
      } else {
        ops_partials.edge3_.partials_ = x_val.transpose() * theta_derivative;
      }
    }
    if (!is_constant_all<T_x_scalar>::value) {
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
        ops_partials.edge2_.partials_[0] = sum(theta_derivative);
      }
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta>
inline return_type_t<T_x, T_beta, T_alpha> bernoulli_logit_glm_lpmf(
    const T_y &y, const T_x &x, const T_alpha &alpha, const T_beta &beta) {
  return bernoulli_logit_glm_lpmf<false>(y, x, alpha, beta);
}
}  // namespace math
}  // namespace stan
#endif
