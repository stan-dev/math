#ifndef STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_BERNOULLI_LOGIT_GLM_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <string>

namespace stan {
  namespace math {

    /**
     * Returns the log PMF of the Generalized Linear Model (GLM)
     * with Bernoulli distribution and logit link function.
     * If containers are supplied, returns the log sum of the probabilities.
     * @tparam T_n type of binary vector of dependent variables (labels);
     * this can also be a single binary value;
     * @tparam T_x type of the matrix of independent variables (features); this
     * should be an Eigen::Matrix type whose number of rows should match the 
     * length of n and whose number of columns should match the length of beta
     * @tparam T_beta type of the weight vector;
     * this can also be a single double value;
     * @tparam T_alpha type of the intercept;
     * this can either be a vector of doubles of a single double value;
     * @param n binary vector parameter
     * @param x design matrix
     * @param beta weight vector
     * @param alpha intercept (in log odds)
     * @return log probability or log sum of probabilities
     * @throw std::domain_error if theta is infinite.
     * @throw std::invalid_argument if container sizes mismatch.
     */
    template <bool propto, typename T_n, typename T_x, typename T_beta,
              typename T_alpha>
    typename return_type<T_x, T_beta, T_alpha>::type
    bernoulli_logit_glm_lpmf(const T_n &n, const T_x &x, const T_beta &beta,
                             const T_alpha &alpha) {
      static const std::string function = "bernoulli_logit_glm_lpmf";
      typedef typename stan::partials_return_type<T_n, T_x, T_beta,
                                                  T_alpha>::type
        T_partials_return;

      using stan::is_constant_struct;
      using std::exp;
      using Eigen::Dynamic;
      using Eigen::Matrix;

      if (!(stan::length(n) && stan::length(x) && stan::length(beta)))
        return 0.0;

      T_partials_return logp(0.0);

      check_bounded(function, "Vector of dependent variables", n, 0, 1);
      check_not_nan(function, "Matrix of independent variables", x);
      check_not_nan(function, "Weight vector", beta);
      check_not_nan(function, "Intercept", alpha);
      check_consistent_sizes(function, "Rows in matrix of independent variables",
                             x.col(0), "Vector of dependent variables",  n);
      check_consistent_sizes(function, "Columns in matrix of independent variables",
                             x.row(0), "Weight vector",  beta);

      if (!include_summand<propto, T_x, T_beta, T_alpha>::value)
        return 0.0;

      scalar_seq_view<T_n> n_vec(n);
      scalar_seq_view<T_beta> beta_vec(beta);

      const size_t N = x.col(0).size();
      const size_t M = x.row(0).size();

      operands_and_partials<T_x, T_beta, T_alpha> ops_partials(x, beta, alpha);

      Matrix<T_partials_return, Dynamic, 1> beta_dbl(M, 1);
      Matrix<T_partials_return, Dynamic, Dynamic> x_dbl(N, M);
      for (size_t m = 0; m < M; ++m) {
        beta_dbl[m] = value_of(beta_vec[m]);
      }

      for (size_t n = 0; n < N; ++n) { // could we vectorise this loop?
        for (size_t m = 0; m < M; ++m) {
          x_dbl(n, m) = value_of(x(n, m));
        }
        const T_partials_return theta_dbl = (x_dbl.row(n) * beta_dbl)[0] +
                                            value_of(alpha);

        const int sign = 2 * n_vec[n] - 1;
        const T_partials_return ntheta = sign * theta_dbl;
        const T_partials_return exp_m_ntheta = exp(-ntheta);

        // Handle extreme values gracefully using Taylor approximations.
        static const double cutoff = 20.0;
        if (ntheta > cutoff)
          logp -= exp_m_ntheta;
        else if (ntheta < -cutoff)
          logp += ntheta;
        else
          logp -= log1p(exp_m_ntheta);

        if (!(is_constant_struct<T_x>::value && is_constant_struct<T_beta>::value
              && is_constant_struct<T_alpha>::value)) {
          T_partials_return theta_derivative;
          if (ntheta > cutoff)
            theta_derivative = -exp_m_ntheta;
          else if (ntheta < -cutoff)
            theta_derivative = sign;
          else
            theta_derivative = sign * exp_m_ntheta / (exp_m_ntheta + 1);
          if (!is_constant_struct<T_beta>::value) {
            ops_partials.edge2_.partials_.col(0).noalias() += theta_derivative *
                                                    x_dbl.row(n).transpose();
          }
          if (!is_constant_struct<T_x>::value) {
            ops_partials.edge1_.partials_.row(n).noalias() += theta_derivative *
                                                              beta_dbl.transpose();
          }
          if (!is_constant_struct<T_alpha>::value) {
            ops_partials.edge3_.partials_[n] += theta_derivative;
          }
        }
      }
      return ops_partials.build(logp);
    }

    template <typename T_n, typename T_x, typename T_beta, typename T_alpha>
    inline
        typename return_type<T_x, T_beta, T_alpha>::type
        bernoulli_logit_glm_lpmf(const T_n &n, const T_x &x, const T_beta &beta,
                                 const T_alpha &alpha) {
      return bernoulli_logit_glm_lpmf<false>(n, x, beta, alpha);
    }
  }
}
#endif
