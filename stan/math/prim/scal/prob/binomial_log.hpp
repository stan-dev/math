#ifndef STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BINOMIAL_LOG_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/log_inv_logit.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/scal/fun/lbeta.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/inc_beta.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
  namespace math {

    // Binomial(n|N, theta)  [N >= 0;  0 <= n <= N;  0 <= theta <= 1]
    template <bool propto,
              typename T_n,
              typename T_N,
              typename T_prob>
    typename return_type<T_prob>::type
    binomial_log(const T_n& n,
                 const T_N& N,
                 const T_prob& theta) {
      typedef typename stan::partials_return_type<T_n, T_N, T_prob>::type
        T_partials_return;

      static const char* function("binomial_log");

      if (!(stan::length(n)
            && stan::length(N)
            && stan::length(theta)))
        return 0.0;

      T_partials_return logp = 0;
      check_bounded(function, "Successes variable", n, 0, N);
      check_nonnegative(function, "Population size parameter", N);
      check_finite(function, "Probability parameter", theta);
      check_bounded(function, "Probability parameter", theta, 0.0, 1.0);
      check_consistent_sizes(function,
                             "Successes variable", n,
                             "Population size parameter", N,
                             "Probability parameter", theta);

      if (!include_summand<propto, T_prob>::value)
        return 0.0;

      VectorView<const T_n> n_vec(n);
      VectorView<const T_N> N_vec(N);
      VectorView<const T_prob> theta_vec(theta);
      size_t size = max_size(n, N, theta);

      OperandsAndPartials<T_prob> operands_and_partials(theta);

      if (include_summand<propto>::value) {
        for (size_t i = 0; i < size; ++i)
          logp += binomial_coefficient_log(N_vec[i], n_vec[i]);
      }

      VectorBuilder<true, T_partials_return, T_prob> log1m_theta(length(theta));
      for (size_t i = 0; i < length(theta); ++i)
        log1m_theta[i] = log1m(value_of(theta_vec[i]));

      for (size_t i = 0; i < size; ++i)
        logp += multiply_log(n_vec[i], value_of(theta_vec[i]))
          + (N_vec[i] - n_vec[i]) * log1m_theta[i];

      if (length(theta) == 1) {
        T_partials_return temp1 = 0;
        T_partials_return temp2 = 0;
        for (size_t i = 0; i < size; ++i) {
          temp1 += n_vec[i];
          temp2 += N_vec[i] - n_vec[i];
        }
        if (!is_constant_struct<T_prob>::value) {
          operands_and_partials.d_x1[0]
            += temp1 / value_of(theta_vec[0])
            - temp2 / (1.0 - value_of(theta_vec[0]));
        }
      } else {
        if (!is_constant_struct<T_prob>::value) {
          for (size_t i = 0; i < size; ++i)
            operands_and_partials.d_x1[i]
              += n_vec[i] / value_of(theta_vec[i])
              - (N_vec[i] - n_vec[i]) / (1.0 - value_of(theta_vec[i]));
        }
      }

      return operands_and_partials.value(logp);
    }

    template <typename T_n,
              typename T_N,
              typename T_prob>
    inline
    typename return_type<T_prob>::type
    binomial_log(const T_n& n,
                 const T_N& N,
                 const T_prob& theta) {
      return binomial_log<false>(n, N, theta);
    }

  }
}
#endif
