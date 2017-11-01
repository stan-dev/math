#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_MIX_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_MIX_HPP

#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/meta/is_constant_struct.hpp>
#include <stan/math/prim/mat/fun/log_sum_exp.hpp>
#include <stan/math/prim/mat/fun/log.hpp>
#include <stan/math/prim/mat/fun/exp.hpp>
#include <stan/math/prim/mat/fun/max.hpp>
#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <string>
#include <algorithm>

namespace stan {
  namespace math {

    /**
     * Return the log mixture density with specified mixing proportion
     * and log densities.
     *
     * \f[
     * \frac{\partial }{\partial p_x}
     *  = \frac{e^{d_x-a}}{e^{d_1-a}p_1 + \cdot\cdot\cdot+e^{d_m-a}p_m}
     * \f]
     *
     * \f[
     * \frac{\partial }{\partial d_x}
     *  = \frac{e^{d_x-a}p_x}{e^{d_1-a}p_1 + \cdot\cdot\cdot+e^{d_m-a}p_m}
     * \f]
     *
     * \f[
     * \mbox{where } a=max(d)
     * \f]
     *
     * @param theta vector of mixing proportions in [0, 1].
     * @param lambda vector of log densities.
     * @return log mixture of densities in specified proportion
     */
    template <typename T_theta, typename T_lam>
    typename return_type<T_theta, T_lam>::type
    log_mix(const T_theta& theta,
            const T_lam& lambda) {
      static const std::string function = "log_mix";
      typedef typename stan::partials_return_type<T_theta, T_lam>::type
        T_partials_return;

      const size_t N = theta.size();

      for (size_t n = 0; n < N; ++n) {
        check_not_nan(function, "lambda", lambda[n]);
        check_not_nan(function, "theta", theta[n]);
        check_finite(function, "lambda", lambda[n]);
        check_finite(function, "theta", theta[n]);
        check_bounded(function, "theta", theta[n], 0, 1);
      }
      check_consistent_sizes(function, "theta", theta, "lambda", lambda);

      Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> theta_dbl(N, 1);
      scalar_seq_view<T_theta> theta_vec(theta);
        for (size_t n = 0; n < N; ++n)
          theta_dbl[n] = value_of(theta_vec[n]);

      Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> lambda_dbl(N, 1);
      scalar_seq_view<T_lam> lam_vec(lambda);
        for (size_t n = 0; n < N; ++n)
          lambda_dbl[n] = value_of(lam_vec[n]);

      /**
       * Vector containing log-prob of each density, log_sum_exp applied as
       * part of return call.
       */
      Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> logp_tmp(N, 1);
        logp_tmp = log(theta_dbl) + lambda_dbl;

      /**
       * Calculate derivatives
       *
       * Exp-normalise to prevent overflow, as: 
       * (exp(x-a) * exp(a)) / (exp(y-a) * exp(a)) = exp(x-a) / exp(y-a)
       */
      Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> exp_lambda_dbl(N, 1);
      double max_val = max(lambda_dbl);
        for (size_t n = 0; n < N; ++n)
          exp_lambda_dbl[n] = exp((lambda_dbl[n] - max_val));

      T_partials_return dot_exp_lam_theta = dot_product(exp_lambda_dbl,
                                                        theta_dbl);

      Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> theta_deriv(N, 1);
          theta_deriv = exp_lambda_dbl / dot_exp_lam_theta;

      Eigen::Matrix<T_partials_return, Eigen::Dynamic, 1> lam_deriv(N, 1);
        for (size_t n = 0; n < N; ++n)
          lam_deriv[n] = theta_deriv[n] * theta_dbl[n];

      operands_and_partials<T_theta, T_lam> ops_partials(theta, lambda);
      if (!(is_constant_struct<T_theta>::value
                                && is_constant_struct<T_lam>::value)) {
        if (!is_constant_struct<T_theta>::value)
          ops_partials.edge1_.partials_ = theta_deriv;

        if (!is_constant_struct<T_lam>::value)
          ops_partials.edge2_.partials_ = lam_deriv;
      }

      return ops_partials.build(log_sum_exp(logp_tmp));
    }
  }
}
#endif
