#ifndef STAN_MATH_PRIM_SCAL_PROB_POISSON_CDF_LOG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_POISSON_CDF_LOG_HPP

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
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <limits>

namespace stan {

  namespace math {

    template <typename T_n, typename T_rate>
    typename return_type<T_rate>::type
    poisson_cdf_log(const T_n& n, const T_rate& lambda) {
      static const char* function("stan::math::poisson_cdf_log");
      typedef typename stan::partials_return_type<T_n, T_rate>::type
        T_partials_return;

      using stan::math::check_not_nan;
      using stan::math::check_nonnegative;
      using stan::math::value_of;
      using stan::math::check_consistent_sizes;

      // Ensure non-zero argument slengths
      if (!(stan::length(n) && stan::length(lambda)))
        return 0.0;

      T_partials_return P(0.0);

      // Validate arguments
      check_not_nan(function, "Rate parameter", lambda);
      check_nonnegative(function, "Rate parameter", lambda);
      check_consistent_sizes(function,
                             "Random variable", n,
                             "Rate parameter", lambda);

      // Wrap arguments into vector views
      VectorView<const T_n> n_vec(n);
      VectorView<const T_rate> lambda_vec(lambda);
      size_t size = max_size(n, lambda);

      // Compute vectorized cdf_log and gradient
      using stan::math::value_of;
      using boost::math::gamma_q;
      using boost::math::lgamma;
      using std::log;
      using std::exp;

      OperandsAndPartials<T_rate> operands_and_partials(lambda);

      // Explicit return for extreme values
      // The gradients are technically ill-defined, but treated as neg infinity
      for (size_t i = 0; i < stan::length(n); i++) {
        if (value_of(n_vec[i]) < 0)
          return operands_and_partials.value(stan::math::negative_infinity());
      }

      for (size_t i = 0; i < size; i++) {
        // Explicit results for extreme values
        // The gradients are technically ill-defined, but treated as zero
        if (value_of(n_vec[i]) == std::numeric_limits<int>::max())
          continue;

        const T_partials_return n_dbl = value_of(n_vec[i]);
        const T_partials_return lambda_dbl = value_of(lambda_vec[i]);
        const T_partials_return log_Pi = log(gamma_q(n_dbl+1, lambda_dbl));

        P += log_Pi;

        if (!is_constant_struct<T_rate>::value)
          operands_and_partials.d_x1[i] += - exp(n_dbl * log(lambda_dbl)
            - lambda_dbl - lgamma(n_dbl+1) - log_Pi);
      }

      return operands_and_partials.value(P);
    }
  }
}
#endif
