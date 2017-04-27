#ifndef STAN_MATH_PRIM_SCAL_PROB_POISSON_LOG_CDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_POISSON_LOG_CDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/gamma_q.hpp>
#include <stan/math/prim/scal/fun/tgamma.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <limits>

namespace stan {
  namespace math {

    // Poisson Log CDF
    template <typename T_n, typename T_rate>
    typename return_type<T_rate>::type
    poisson_log_cdf(const T_n& n, const T_rate& log_lambda) {
      static const char* function("poisson_log_cdf");
      typedef typename stan::partials_return_type<T_n, T_rate>::type
        T_partials_return;

      if (!(stan::length(n) && stan::length(log_lambda)))
        return 1.0;

      T_partials_return P(1.0);

      check_not_nan(function, "Log rate parameter", log_lambda);
      check_consistent_sizes(function,
                             "Random variable", n,
                             "Log rate parameter", log_lambda);

      VectorView<const T_n> n_vec(n);
      VectorView<const T_rate> log_lambda_vec(log_lambda);
      size_t size = max_size(n, log_lambda);

      using std::exp;
      using std::pow;

      OperandsAndPartials<T_rate> operands_and_partials(log_lambda);

      // Explicit return for extreme values
      // The gradients are technically ill-defined, but treated as zero
      for (size_t i = 0; i < stan::length(n); i++) {
        if (value_of(n_vec[i]) < 0)
          return operands_and_partials.value(0.0);
      }

      for (size_t i = 0; i < size; i++) {
        // Explicit results for extreme values
        // The gradients are technically ill-defined, but treated as zero
        if (value_of(n_vec[i]) == std::numeric_limits<int>::max())
          continue;

        const T_partials_return n_dbl = value_of(n_vec[i]);
        const T_partials_return log_lambda_dbl = value_of(log_lambda_vec[i]);
        const T_partials_return lambda_dbl = exp(log_lambda_dbl);
        const T_partials_return Pi = gamma_q(n_dbl+1, lambda_dbl);

        P *= Pi;

        if (!is_constant_struct<T_rate>::value)
          operands_and_partials.d_x1[i] -= exp(-lambda_dbl)
            * pow(lambda_dbl, n_dbl) / tgamma(n_dbl+1) / Pi;
      }

      if (!is_constant_struct<T_rate>::value) {
        for (size_t i = 0; i < stan::length(log_lambda); ++i) {
        const T_partials_return log_lambda_dbl = value_of(log_lambda_vec[i]);
        const T_partials_return lambda_dbl = exp(log_lambda_dbl);
          operands_and_partials.d_x1[i] *= P * lambda_dbl;
        }
      }
      return operands_and_partials.value(P);
    }

  }
}
#endif
