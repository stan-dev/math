#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_MIX_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_MIX_HPP

#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/fun/log.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/is_constant_struct.hpp>

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
     * \f[
     * \frac{\partial}{\partial \theta}
     * \mbox{log\_mix}(\theta, \lambda_1, \lambda_2)
     * = FIXME
     * \f]
     *
     * \f[
     * \frac{\partial}{\partial \lambda_1}
     * \mbox{log\_mix}(\theta, \lambda_1, \lambda_2)
     * = FIXME
     * \f]
     *
     * \f[
     * \frac{\partial}{\partial \lambda_2}
     * \mbox{log\_mix}(\theta, \lambda_1, \lambda_2)
     * = FIXME
     * \f]
     *
     * @param[in] theta mixing proportion in [0, 1].
     * @param lambda1 first log density.
     * @param lambda2 second log density.
     * @return log mixture of densities in specified proportion
     */
    template <typename T_theta, typename T_lam>
    typename return_type<T_theta, T_lam>::type
    log_mix(const T_theta& theta,
            const T_lam& lambda) {
    typedef typename stan::partials_return_type<T_theta, T_lam>::type
      T_partials_return;

    const size_t N = theta.size();

    scalar_seq_view<T_theta> theta_vec(theta);
    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> theta_dbl(N, 1);
      for (size_t n = 0; n < N; ++n) {
        theta_dbl[n] = value_of(theta_vec[n]);
      }

    scalar_seq_view<T_lam> lam_vec(lambda);
    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> lambda_dbl(N, 1);
      for (size_t n = 0; n < N; ++n) {
        lambda_dbl[n] = value_of(lam_vec[n]);
      }

    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> exp_lambda_dbl(N, 1);
        exp_lambda_dbl = exp(lambda_dbl);

    T_partials_return logp(0.0);
    T_partials_return logp_tmp(0.0);

    for (size_t n = 0; n < N; ++n) 
      logp_tmp += theta_dbl[n] * exp_lambda_dbl[n];

    logp += log(logp_tmp);

    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> theta_deriv(N, 1);
        theta_deriv = exp_lambda_dbl / dot_product(exp_lambda_dbl, theta_dbl); 

    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> lam_deriv(N, 1);
      for (size_t n = 0; n < N; ++n) 
        lam_deriv[n] = theta_deriv[n] * theta_dbl[n];


      operands_and_partials<T_theta, T_lam> ops_partials(theta, lambda);
      if (!(is_constant_struct<T_theta>::value && is_constant_struct<T_lam>::value)) {
        if (!is_constant_struct<T_theta>::value) {
          ops_partials.edge1_.partials_ = theta_deriv;
        }
        if (!is_constant_struct<T_lam>::value) {
          ops_partials.edge2_.partials_ = lam_deriv;
        }
      }
      return ops_partials.build(logp);
    }
  }
}
#endif
