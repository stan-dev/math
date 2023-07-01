#ifndef STAN_MATH_PRIM_PROB_NORMAL_ID_GLM_LPDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_ID_GLM_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/isfinite.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Returns the log PDF of the Generalized Linear Model (GLM)
 * with Normal distribution and id link function.
 * If containers are supplied, returns the log sum of the probabilities.
 * The idea is that normal_id_glm_lpdf(y, x, alpha, beta, sigma) should
 * compute a more efficient version of normal_lpdf(y, alpha + x * beta, sigma)
 * by using analytically simplified gradients.
 *
 * @tparam T_y type of vector of dependent variables (labels);
 * @tparam T_x type of the matrix of independent variables (features)
 * @tparam T_alpha type of the intercept(s);
 * this can be a vector (of the same length as y) of intercepts or a single
 * value (for models with constant intercept);
 * @tparam T_beta type of the weight vector;
 * this can also be a single value;
 * @tparam T_scale type of the (positive) scale(s);
 * this can be a vector (of the same length as y, for heteroskedasticity)
 * or a scalar.
 *
 * @param y scalar or vector of dependent variables. If it is a scalar it will
 * be broadcast - used for all instances.
 * @param x design matrix or row vector. If it is a row vector it will be
 * broadcast - used for all instances.
 * @param alpha intercept (in log odds)
 * @param beta weight vector
 * @param sigma (Sequence of) scale parameters for the normal
 * distribution.
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if x, beta or alpha is infinite.
 * @throw std::domain_error if the scale is not positive.
 * @throw std::invalid_argument if container sizes mismatch.
 */
template <bool propto, typename T_y, typename T_x, typename T_alpha,
          typename T_beta, typename T_scale, require_matrix_t<T_x>* = nullptr>
return_type_t<T_y, T_x, T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta,
    const T_scale& sigma) {
  using Eigen::Array;
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using std::isfinite;
  constexpr int T_x_rows = T_x::RowsAtCompileTime;
  using T_partials_return
      = partials_return_t<T_y, T_x, T_alpha, T_beta, T_scale>;
  using T_scale_val = typename std::conditional_t<
      is_vector<T_scale>::value,
      Eigen::Array<partials_return_t<T_scale>, -1, 1>,
      partials_return_t<T_scale>>;
  using T_y_scaled_tmp =
      typename std::conditional_t<T_x_rows == 1, T_partials_return,
                                  Array<T_partials_return, Dynamic, 1>>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_x_ref = ref_type_if_t<!is_constant<T_x>::value, T_x>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_alpha>::value, T_alpha>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_beta>::value, T_beta>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;

  const size_t N_instances = T_x_rows == 1 ? stan::math::size(y) : x.rows();
  const size_t N_attributes = x.cols();

  static const char* function = "normal_id_glm_lpdf";
  check_consistent_size(function, "Vector of dependent variables", y,
                        N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_consistent_size(function, "Vector of scale parameters", sigma,
                        N_instances);
  check_consistent_size(function, "Vector of intercepts", alpha, N_instances);
  T_sigma_ref sigma_ref = sigma;
  const auto& sigma_val = value_of(sigma_ref);
  const auto& sigma_val_vec = to_ref(as_column_vector_or_scalar(sigma_val));
  check_positive_finite(function, "Scale vector", sigma_val_vec);

  if (size_zero(y, sigma)) {
    return 0;
  }
  if (!include_summand<propto, T_y, T_x, T_alpha, T_beta, T_scale>::value) {
    return 0;
  }

  T_y_ref y_ref = y;
  T_x_ref x_ref = x;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  const auto& y_val = value_of(y_ref);
  const auto& x_val = to_ref_if<!is_constant<T_beta>::value>(value_of(x_ref));
  const auto& alpha_val = value_of(alpha_ref);
  const auto& beta_val = value_of(beta_ref);

  const auto& y_val_vec = as_column_vector_or_scalar(y_val);
  const auto& alpha_val_vec = as_column_vector_or_scalar(alpha_val);
  const auto& beta_val_vec = to_ref_if<!is_constant<T_x>::value>(
      as_column_vector_or_scalar(beta_val));

  T_scale_val inv_sigma = 1.0 / as_array_or_scalar(sigma_val_vec);

  // the most efficient way to calculate this depends on template parameters
  T_partials_return y_scaled_sq_sum;

  Array<T_partials_return, Dynamic, 1> y_scaled(N_instances);
  if (T_x_rows == 1) {
    T_y_scaled_tmp y_scaled_tmp
        = forward_as<T_y_scaled_tmp>((x_val * beta_val_vec).coeff(0, 0));
    y_scaled = (as_array_or_scalar(y_val_vec) - y_scaled_tmp
                - as_array_or_scalar(alpha_val_vec))
               * inv_sigma;
  } else {
    y_scaled = x_val * beta_val_vec;
    y_scaled = (as_array_or_scalar(y_val_vec) - y_scaled
                - as_array_or_scalar(alpha_val_vec))
               * inv_sigma;
  }

  auto ops_partials
      = make_partials_propagator(y_ref, x_ref, alpha_ref, beta_ref, sigma_ref);

  if (!(is_constant_all<T_y, T_x, T_beta, T_alpha, T_scale>::value)) {
    Matrix<T_partials_return, Dynamic, 1> mu_derivative = inv_sigma * y_scaled;
    if (!is_constant_all<T_y>::value) {
      if (is_vector<T_y>::value) {
        partials<0>(ops_partials) = -mu_derivative;
      } else {
        partials<0>(ops_partials)[0] = -mu_derivative.sum();
      }
    }
    if (!is_constant_all<T_x>::value) {
      if (T_x_rows == 1) {
        edge<1>(ops_partials).partials_
            = forward_as<Array<T_partials_return, Dynamic, T_x_rows>>(
                beta_val_vec * sum(mu_derivative));
      } else {
        edge<1>(ops_partials).partials_
            = (beta_val_vec * mu_derivative.transpose()).transpose();
      }
    }
    if (!is_constant_all<T_beta>::value) {
      if (T_x_rows == 1) {
        edge<3>(ops_partials).partials_
            = forward_as<Matrix<T_partials_return, 1, Dynamic>>(
                mu_derivative.sum() * x_val);
      } else {
        partials<3>(ops_partials) = mu_derivative.transpose() * x_val;
      }
    }
    if (!is_constant_all<T_alpha>::value) {
      if (is_vector<T_alpha>::value) {
        partials<2>(ops_partials) = mu_derivative;
      } else {
        partials<2>(ops_partials)[0] = sum(mu_derivative);
      }
    }
    if (!is_constant_all<T_scale>::value) {
      if (is_vector<T_scale>::value) {
        Array<T_partials_return, Dynamic, 1> y_scaled_sq = y_scaled * y_scaled;
        y_scaled_sq_sum = sum(y_scaled_sq);
        partials<4>(ops_partials) = (y_scaled_sq - 1) * inv_sigma;
      } else {
        y_scaled_sq_sum = sum(y_scaled * y_scaled);
        partials<4>(ops_partials)[0]
            = (y_scaled_sq_sum - N_instances)
              * forward_as<partials_return_t<T_sigma_ref>>(inv_sigma);
      }
    } else {
      y_scaled_sq_sum = sum(y_scaled * y_scaled);
    }
  } else {
    y_scaled_sq_sum = sum(y_scaled * y_scaled);
  }

  if (!isfinite(y_scaled_sq_sum)) {
    check_finite(function, "Vector of dependent variables", y_val_vec);
    check_finite(function, "Weight vector", beta_val_vec);
    check_finite(function, "Intercept", alpha_val_vec);
    // if all other checks passed, next will only fail if x is not finite
    check_finite(function, "Matrix of independent variables", y_scaled_sq_sum);
  }

  // Compute log probability.
  T_partials_return logp(0.0);
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N_instances;
  }
  if (include_summand<propto, T_scale>::value) {
    if (is_vector<T_scale>::value) {
      logp -= sum(log(sigma_val_vec));
    } else {
      logp -= N_instances
              * log(forward_as<partials_return_t<T_sigma_ref>>(sigma_val_vec));
    }
  }
  logp -= 0.5 * y_scaled_sq_sum;

  return ops_partials.build(logp);
}

template <typename T_y, typename T_x, typename T_alpha, typename T_beta,
          typename T_scale>
inline return_type_t<T_y, T_x, T_alpha, T_beta, T_scale> normal_id_glm_lpdf(
    const T_y& y, const T_x& x, const T_alpha& alpha, const T_beta& beta,
    const T_scale& sigma) {
  return normal_id_glm_lpdf<false>(y, x, alpha, beta, sigma);
}
}  // namespace math
}  // namespace stan
#endif
