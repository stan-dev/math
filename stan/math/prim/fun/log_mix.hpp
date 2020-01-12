#ifndef STAN_MATH_PRIM_FUN_LOG_MIX_HPP
#define STAN_MATH_PRIM_FUN_LOG_MIX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/value_of.hpp>
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
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(double theta, double lambda1, double lambda2) {
  using std::log;
  check_not_nan("log_mix", "lambda1", lambda1);
  check_not_nan("log_mix", "lambda2", lambda2);
  check_bounded("log_mix", "theta", theta, 0, 1);
  return log_sum_exp(log(theta) + lambda1, log1m(theta) + lambda2);
}

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(double theta, double lambda1, int lambda2) {
  return log_mix(theta, lambda1, static_cast<double>(lambda2));
}

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(double theta, int lambda1, double lambda2) {
  return log_mix(theta, static_cast<double>(lambda1), lambda2);
}

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(int theta, double lambda1, double lambda2) {
  return log_mix(static_cast<double>(theta), lambda1, lambda2);
}

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(double theta, int lambda1, int lambda2) {
  return log_mix(theta, static_cast<double>(lambda1),
                 static_cast<double>(lambda2));
}

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(int theta, double lambda1, int lambda2) {
  return log_mix(static_cast<double>(theta), lambda1,
                 static_cast<double>(lambda2));
}

/**
 * Return the log mixture density with specified mixing proportion
 * and log densities.
 *
 * @param[in] theta mixing proportion in [0, 1].
 * @param[in] lambda1 first log density.
 * @param[in] lambda2 second log density.
 * @return log mixture of densities in specified proportion
 */
inline double log_mix(int theta, int lambda1, int lambda2) {
  return log_mix(static_cast<double>(theta), static_cast<double>(lambda1),
                 static_cast<double>(lambda2));
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
 * @param theta vector of mixing proportions in [0, 1].
 * @param lambda vector of log densities.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta, typename T_lam>
return_type_t<T_theta, T_lam> log_mix(const T_theta& theta,
                                      const T_lam& lambda) {
  static const char* function = "log_mix";
  using T_partials_return = partials_return_t<T_theta, T_lam>;

  using T_partials_vec = typename Eigen::Matrix<T_partials_return, -1, 1>;

  const int N = size(theta);

  check_bounded(function, "theta", theta, 0, 1);
  check_not_nan(function, "lambda", lambda);
  check_not_nan(function, "theta", theta);
  check_finite(function, "lambda", lambda);
  check_finite(function, "theta", theta);
  check_consistent_sizes(function, "theta", theta, "lambda", lambda);

  scalar_seq_view<T_theta> theta_vec(theta);
  T_partials_vec theta_dbl(N);
  for (int n = 0; n < N; ++n) {
    theta_dbl[n] = value_of(theta_vec[n]);
  }

  scalar_seq_view<T_lam> lam_vec(lambda);
  T_partials_vec lam_dbl(N);
  for (int n = 0; n < N; ++n) {
    lam_dbl[n] = value_of(lam_vec[n]);
  }

  T_partials_return logp = log_sum_exp((log(theta_dbl) + lam_dbl).eval());

  T_partials_vec theta_deriv(N);
  theta_deriv.array() = (lam_dbl.array() - logp).exp();

  T_partials_vec lam_deriv = theta_deriv.cwiseProduct(theta_dbl);

  operands_and_partials<T_theta, T_lam> ops_partials(theta, lambda);
  if (!is_constant_all<T_theta>::value) {
    for (int n = 0; n < N; ++n) {
      ops_partials.edge1_.partials_[n] = theta_deriv[n];
    }
  }

  if (!is_constant_all<T_lam>::value) {
    for (int n = 0; n < N; ++n) {
      ops_partials.edge2_.partials_[n] = lam_deriv[n];
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
 * @param theta vector of mixing proportions in [0, 1].
 * @param lambda array containing vectors of log densities.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta, typename T_lam, int R, int C>
return_type_t<T_theta, std::vector<Eigen::Matrix<T_lam, R, C>>> log_mix(
    const T_theta& theta,
    const std::vector<Eigen::Matrix<T_lam, R, C>>& lambda) {
  static const char* function = "log_mix";
  using T_partials_return
      = partials_return_t<T_theta, std::vector<Eigen::Matrix<T_lam, R, C>>>;

  using T_partials_vec = typename Eigen::Matrix<T_partials_return, -1, 1>;

  using T_partials_mat = typename Eigen::Matrix<T_partials_return, -1, -1>;

  using T_lamvec_type = typename std::vector<Eigen::Matrix<T_lam, R, C>>;

  const int N = size(lambda);
  const int M = theta.size();

  check_bounded(function, "theta", theta, 0, 1);
  check_not_nan(function, "theta", theta);
  check_finite(function, "theta", theta);
  for (int n = 0; n < N; ++n) {
    check_not_nan(function, "lambda", lambda[n]);
    check_finite(function, "lambda", lambda[n]);
    check_consistent_sizes(function, "theta", theta, "lambda", lambda[n]);
  }

  scalar_seq_view<T_theta> theta_vec(theta);
  T_partials_vec theta_dbl(M);
  for (int m = 0; m < M; ++m) {
    theta_dbl[m] = value_of(theta_vec[m]);
  }

  T_partials_mat lam_dbl(M, N);
  vector_seq_view<T_lamvec_type> lam_vec(lambda);
  for (int n = 0; n < N; ++n) {
    for (int m = 0; m < M; ++m) {
      lam_dbl(m, n) = value_of(lam_vec[n][m]);
    }
  }

  T_partials_mat logp_tmp = log(theta_dbl).replicate(1, N) + lam_dbl;

  T_partials_vec logp(N);
  for (int n = 0; n < N; ++n) {
    logp[n] = log_sum_exp(logp_tmp.col(n).eval());
  }

  operands_and_partials<T_theta, T_lamvec_type> ops_partials(theta, lambda);

  if (!is_constant_all<T_theta, T_lam>::value) {
    T_partials_mat derivs
        = (lam_dbl - logp.transpose().replicate(M, 1))
              .unaryExpr([](T_partials_return x) { return exp(x); });
    if (!is_constant_all<T_theta>::value) {
      for (int m = 0; m < M; ++m) {
        ops_partials.edge1_.partials_[m] = derivs.row(m).sum();
      }
    }

    if (!is_constant_all<T_lam>::value) {
      for (int n = 0; n < N; ++n) {
        ops_partials.edge2_.partials_vec_[n]
            = derivs.col(n).cwiseProduct(theta_dbl);
      }
    }
  }
  return ops_partials.build(logp.sum());
}
/**
 * Return the log mixture density given specified mixing proportions
 * and array of log density arrays.
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
 * @param theta vector of mixing proportions in [0, 1].
 * @param lambda array containing arrays of log densities.
 * @return log mixture of densities in specified proportion
 */
template <typename T_theta, typename T_lam>
return_type_t<T_theta, std::vector<std::vector<T_lam>>> log_mix(
    const T_theta& theta, const std::vector<std::vector<T_lam>>& lambda) {
  static const char* function = "log_mix";
  using T_partials_return
      = partials_return_t<T_theta, std::vector<std::vector<T_lam>>>;
  using T_partials_vec = typename Eigen::Matrix<T_partials_return, -1, 1>;
  using T_partials_mat = typename Eigen::Matrix<T_partials_return, -1, -1>;
  using T_lamvec_type = typename std::vector<std::vector<T_lam>>;

  const int N = size(lambda);
  const int M = theta.size();

  check_bounded(function, "theta", theta, 0, 1);
  check_not_nan(function, "theta", theta);
  check_finite(function, "theta", theta);
  for (int n = 0; n < N; ++n) {
    check_not_nan(function, "lambda", lambda[n]);
    check_finite(function, "lambda", lambda[n]);
    check_consistent_sizes(function, "theta", theta, "lambda", lambda[n]);
  }

  scalar_seq_view<T_theta> theta_vec(theta);
  T_partials_vec theta_dbl(M);
  for (int m = 0; m < M; ++m) {
    theta_dbl[m] = value_of(theta_vec[m]);
  }

  T_partials_mat lam_dbl(M, N);
  for (int n = 0; n < N; ++n) {
    for (int m = 0; m < M; ++m) {
      lam_dbl(m, n) = value_of(lambda[n][m]);
    }
  }

  T_partials_mat logp_tmp = log(theta_dbl).replicate(1, N) + lam_dbl;

  T_partials_vec logp(N);
  for (int n = 0; n < N; ++n) {
    logp[n] = log_sum_exp(logp_tmp.col(n).eval());
  }

  T_partials_mat derivs
      = (lam_dbl - logp.transpose().replicate(M, 1))
            .unaryExpr([](T_partials_return x) { return exp(x); });

  T_partials_mat lam_deriv(M, N);
  for (int n = 0; n < N; ++n) {
    lam_deriv.col(n) = derivs.col(n).cwiseProduct(theta_dbl);
  }

  operands_and_partials<T_theta, T_lamvec_type> ops_partials(theta, lambda);
  if (!is_constant_all<T_theta>::value) {
    for (int m = 0; m < M; ++m) {
      ops_partials.edge1_.partials_[m] = derivs.row(m).sum();
    }
  }

  if (!is_constant_all<T_lam>::value) {
    for (int n = 0; n < N; ++n) {
      for (int m = 0; m < M; ++m) {
        ops_partials.edge2_.partials_vec_[n][m] = lam_deriv(m, n);
      }
    }
  }
  return ops_partials.build(logp.sum());
}
}  // namespace math
}  // namespace stan
#endif
