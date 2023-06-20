#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
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
 *
 * @tparam T_y type of binary vector of dependent variables (labels);
 * this can also be a single binary value;
 * @tparam T_x type of the matrix of independent variables (features)
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector
 *
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
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, require_matrix_t<T_x>* = nullptr>
return_type_t<T_x, T_alpha, T_beta> bernoulli_logit_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::log1p;
  using Eigen::Matrix;
  using std::exp;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_partials_return = partials_return_t<T_y, T_x, T_alpha, T_beta>;
  using T_ytheta_tmp =
      typename std::conditional_t<T_x_rows == 1, T_partials_return,
                                  Array<T_partials_return, Dynamic, 1>>;
  using T_x_ref = ref_type_if_t<!is_constant<T_x>::value, T_x>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_alpha>::value, T_alpha>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();

  static const char* function = "bernoulli_logit_glm_lpmf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);
  if (size_zero(y)) {
    return 0;
  }

  const auto& y_ref = to_ref(y);
  check_bounded(function, "Vector of dependent variables", y_ref, 0, 1);

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

  const auto& y_val_vec = as_column_vector_or_scalar(y_val);
  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val);
  const auto& beta_val_vec = to_ref_if<!is_constant<T_x>::value>(
      as_column_vector_or_scalar(beta_val));

  auto signs = to_ref_if<!is_constant_all<T_beta, T_x, T_alpha>::value>(
      2 * as_array_or_scalar(y_val_vec) - 1);

  Array<T_partials_return, Dynamic, 1> ytheta(N_instances);
  if (T_x_rows == 1) {
    T_ytheta_tmp ytheta_tmp
        = forward_as<T_ytheta_tmp>((x_val * beta_val_vec)(0, 0));
    ytheta = signs * (ytheta_tmp + as_array_or_scalar(alpha_val_vec));
  } else {
    ytheta = (x_val * beta_val_vec).array();
    ytheta = signs * (ytheta + as_array_or_scalar(alpha_val_vec));
  }

  // Compute the log-density and handle extreme values gracefully
  // using Taylor approximations.
  // And compute the derivatives wrt theta.
  static const double cutoff = 20.0;
  Eigen::Array<T_partials_return, Dynamic, 1> exp_m_ytheta = exp(-ytheta);
  T_partials_return logp = sum(
      (ytheta > cutoff)
          .select(-exp_m_ytheta,
                  (ytheta < -cutoff).select(ytheta, -log1p(exp_m_ytheta))));

  if (!std::isfinite(logp)) {
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Intercept", alpha);
    check_finite(function, "Matrix of independent variables", ytheta);
  }

  auto ops_partials = make_partials_propagator(x_ref, alpha_ref, beta_ref);
  // Compute the necessary derivatives.
  if (!is_constant_all<T_beta, T_x, T_alpha>::value) {
    Matrix<T_partials_return, Dynamic, 1> theta_derivative
        = (ytheta > cutoff)
              .select(-exp_m_ytheta,
                      (ytheta < -cutoff)
                          .select(signs,
                                  signs * exp_m_ytheta / (exp_m_ytheta + 1)));
    if (!is_constant_all<T_beta>::value) {
      if (T_x_rows == 1) {
        edge<2>(ops_partials).partials_
            = forward_as<Matrix<T_partials_return, 1, Dynamic>>(
                theta_derivative.sum() * x_val);
      } else {
        partials<2>(ops_partials) = x_val.transpose() * theta_derivative;
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
      partials<1>(ops_partials) = theta_derivative;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta>
inline return_type_t<T_x, T_beta, T_alpha> bernoulli_logit_glm_lpmf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta) {
  return bernoulli_logit_glm_lpmf<false>(y, x, alpha, beta);
}
}  // namespace math
}  // namespace stan
#endif
