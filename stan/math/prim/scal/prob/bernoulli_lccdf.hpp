#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LCCDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
  namespace math {

    template <typename T_n, typename T_prob>
    typename return_type<T_prob>::type
    bernoulli_lccdf(const T_n& n, const T_prob& theta) {
      static const char* function("bernoulli_lccdf");
      typedef typename stan::partials_return_type<T_n, T_prob>::type
        T_partials_return;

      if (!(stan::length(n) && stan::length(theta)))
        return 0.0;

      T_partials_return P(0.0);

      check_finite(function, "Probability parameter", theta);
      check_bounded(function, "Probability parameter", theta, 0.0, 1.0);
      check_consistent_sizes(function,
                             "Random variable", n,
                             "Probability parameter", theta);

      VectorView<const T_n> n_vec(n);
      VectorView<const T_prob> theta_vec(theta);
      size_t size = max_size(n, theta);

      using std::log;
      OperandsAndPartials<T_prob> operands_and_partials(theta);

      // Explicit return for extreme values
      // The gradients are technically ill-defined, but treated as zero
      for (size_t i = 0; i < stan::length(n); i++) {
        if (value_of(n_vec[i]) < 0)
          return operands_and_partials.value(0.0);
      }

      for (size_t i = 0; i < size; i++) {
        // Explicit results for extreme values
        // The gradients are technically ill-defined, but treated as zero
        if (value_of(n_vec[i]) >= 1) {
          return operands_and_partials.value(negative_infinity());
        } else {
          const T_partials_return Pi = value_of(theta_vec[i]);

          P += log(Pi);

          if (!is_constant_struct<T_prob>::value)
            operands_and_partials.d_x1[i] += 1 / Pi;
        }
      }

      return operands_and_partials.value(P);
    }
  }
}
#endif
