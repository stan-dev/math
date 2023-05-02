#ifndef STAN_MATH_PRIM_FUN_LOG_MIX_HPP
#define STAN_MATH_PRIM_FUN_LOG_MIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/partials_propagator.hpp>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * \f[
 * \mbox{log\_mix}(\theta, \lambda_1, \lambda_2)
 * = \log \left( \theta \lambda_1 + (1 - \theta) \lambda_2 \right).
 * \f]
 *
 * @tparam T_theta type of mixing proportion - must be an arithmetic type
 * @tparam T_lambda1 type of first log density - must be an arithmetic type
 * @tparam T_lambda2 type of second log density - must be an arithmetic type
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta, typename T_lambda1, typename T_lambda2,
          require_all_arithmetic_t<T_theta, T_lambda1, T_lambda2>* = nullptr>
inline double log_mix(T_theta theta, T_lambda1 lambda1, T_lambda2 lambda2) {
  using std::log;
  check_not_nan("log_mix", "lambda1", lambda1);
  check_not_nan("log_mix", "lambda2", lambda2);
  check_bounded("log_mix", "theta", theta, 0, 1);
  return log_sum_exp(log(theta) + lambda1, log1m(theta) + lambda2);
}

/**
 * Return the log mixture density with specified mixing proportions
 * and log densities.
 *
 * \f[
 * \frac{\partial }{\partial p_x}
 * \log\left(\exp^{\log\left(p_1\right)+d_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+d_n}\right)
 * =\frac{e^{d_x}}{e^{d_1}p_1+\cdot\cdot\cdot+e^{d_m}p_m}
 * \f]
 *
 * \f[
 * \frac{\partial }{\partial d_x}
 * \log\left(\exp^{\log\left(p_1\right)+d_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+d_n}\right)
 * =\frac{e^{d_x}p_x}{e^{d_1}p_1+\cdot\cdot\cdot+e^{d_m}p_m}
 * \f]
 *
 * @tparam T_theta Type of theta. This can be a scalar, std vector or row/column
 * vector.
 * @tparam T_lam Type of lambda. This can be a scalar, std vector or row/column
 * vector.
 * @param theta std/row/col vector of mixing proportions in [0, 1].
 * @param lambda std/row/col vector of log densities.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta, typename T_lam,
          require_any_vector_t<T_theta, T_lam>* = nullptr>
return_type_t<T_theta, T_lam> log_mix(const T_theta& theta,
                                      const T_lam& lambda) {
  static const char* function = "log_mix";
  using T_partials_return = partials_return_t<T_theta, T_lam>;
  using T_partials_vec =
      typename Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using T_theta_ref = ref_type_t<T_theta>;
  using T_lam_ref = ref_type_t<T_lam>;

  check_consistent_sizes(function, "theta", theta, "lambda", lambda);
  T_theta_ref theta_ref = theta;
  T_lam_ref lambda_ref = lambda;
  check_bounded(function, "theta", theta_ref, 0, 1);
  check_finite(function, "lambda", lambda_ref);

  const auto& theta_dbl
      = to_ref(value_of(as_column_vector_or_scalar(theta_ref)));
  const auto& lam_dbl
      = to_ref(value_of(as_column_vector_or_scalar(lambda_ref)));

  T_partials_return logp = log_sum_exp(log(theta_dbl) + lam_dbl);

  auto ops_partials = make_partials_propagator(theta_ref, lambda_ref);
  if (!is_constant_all<T_lam, T_theta>::value) {
    T_partials_vec theta_deriv = (lam_dbl.array() - logp).exp();
    if (!is_constant_all<T_lam>::value) {
      partials<1>(ops_partials) = theta_deriv.cwiseProduct(theta_dbl);
    }
    if (!is_constant_all<T_theta>::value) {
      partials<0>(ops_partials) = std::move(theta_deriv);
    }
  }
  return ops_partials.build(logp);
}

/**
 * Return the log mixture density given specified mixing proportions
 * and array of log density vectors.
 *
 * \f[
 * \frac{\partial }{\partial p_x}\left[
 * \log\left(\exp^{\log\left(p_1\right)+d_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+d_n}\right)+
 * \log\left(\exp^{\log\left(p_1\right)+f_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+f_n}\right)\right]
 * =\frac{e^{d_x}}{e^{d_1}p_1+\cdot\cdot\cdot+e^{d_m}p_m}+
 * \frac{e^{f_x}}{e^{f_1}p_1+\cdot\cdot\cdot+e^{f_m}p_m}
 * \f]
 *
 * \f[
 * \frac{\partial }{\partial d_x}\left[
 * \log\left(\exp^{\log\left(p_1\right)+d_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+d_n}\right)
 * +\log\left(\exp^{\log\left(p_1\right)+f_1}+\cdot\cdot\cdot+
 * \exp^{\log\left(p_n\right)+f_n}\right)\right]
 * =\frac{e^{d_x}p_x}{e^{d_1}p_1+\cdot\cdot\cdot+e^{d_m}p_m}
 * \f]
 *
 * @tparam T_theta Type of theta. This can be a scalar, std vector or row/column
 * vector
 * @tparam T_lam Type of vector in std vector lambda. This can be std vector or
 * row/column vector.
 * @param theta std/row/col vector of mixing proportions in [0, 1].
 * @param lambda std vector containing std/row/col vectors of log densities.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta, typename T_lam, require_vector_t<T_lam>* = nullptr>
return_type_t<T_theta, std::vector<T_lam>> log_mix(
    const T_theta& theta, const std::vector<T_lam>& lambda) {
  static const char* function = "log_mix";
  using T_partials_return = partials_return_t<T_theta, std::vector<T_lam>>;
  using T_partials_vec =
      typename Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1>;
  using T_partials_mat =
      typename Eigen::Matrix<T_partials_return, Eigen::Dynamic, Eigen::Dynamic>;
  using T_theta_ref = ref_type_t<T_theta>;

  const int N = stan::math::size(lambda);
  const int M = theta.size();

  T_theta_ref theta_ref = theta;
  check_bounded(function, "theta", theta_ref, 0, 1);
  for (int n = 0; n < N; ++n) {
    check_not_nan(function, "lambda", lambda[n]);
    check_finite(function, "lambda", lambda[n]);
    check_consistent_sizes(function, "theta", theta, "lambda", lambda[n]);
  }

  const auto& theta_dbl
      = to_ref(value_of(as_column_vector_or_scalar(theta_ref)));

  T_partials_mat lam_dbl(M, N);
  for (int n = 0; n < N; ++n) {
    lam_dbl.col(n) = value_of(as_column_vector_or_scalar(lambda[n]));
  }

  T_partials_mat logp_tmp = lam_dbl.colwise() + log(theta_dbl);
  T_partials_vec logp(N);
  for (int n = 0; n < N; ++n) {
    logp[n] = log_sum_exp(logp_tmp.col(n));
  }

  auto ops_partials = make_partials_propagator(theta_ref, lambda);
  if (!is_constant_all<T_theta, T_lam>::value) {
    T_partials_mat derivs = exp(lam_dbl.rowwise() - logp.transpose());
    if (!is_constant_all<T_theta>::value) {
      partials<0>(ops_partials) = derivs.rowwise().sum();
    }
    if (!is_constant_all<T_lam>::value) {
      for (int n = 0; n < N; ++n) {
        as_column_vector_or_scalar(partials_vec<1>(ops_partials)[n])
            = derivs.col(n).cwiseProduct(theta_dbl);
      }
    }
  }
  return ops_partials.build(logp.sum());
}

}  // namespace math
}  // namespace stan
#endif
