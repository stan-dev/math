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

    T_partials_return logp(0.0);

    const size_t N = theta.size();

    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> theta_dbl(N, 1);
    {
      scalar_seq_view<T_theta> theta_vec(theta);
      for (size_t n = 0; n < N; ++n) {
        theta_dbl[n] = value_of(theta_vec[n]);
      }
    }

    Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> lambda_dbl(N, 1);
    {
      scalar_seq_view<T_lam> lam_vec(lambda);
      for (size_t n = 0; n < N; ++n) {
        lambda_dbl[n] = value_of(lam_vec[n]);
      }
    }
    Eigen::Matrix<typename return_type<T_theta, T_lam>::type, Eigen::Dynamic, 1> lgp(N,1);

    lgp = log(theta) + lambda;

    logp += value_of(log_sum_exp(lgp));

    return logp;
    }

  }
}
#endif
