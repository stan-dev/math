#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_LPDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_LPDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/contains_nonconstant_struct.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
  namespace math {

    // Pareto(y|y_m, alpha)  [y > y_m;  y_m > 0;  alpha > 0]
    template <bool propto,
              typename T_y, typename T_scale, typename T_shape>
    typename return_type<T_y, T_scale, T_shape>::type
    pareto_lpdf(const T_y& y, const T_scale& y_min, const T_shape& alpha) {
      static const char* function("pareto_lpdf");
      typedef typename stan::partials_return_type<T_y, T_scale, T_shape>::type
        T_partials_return;

      using std::log;

      if (!(stan::length(y)
            && stan::length(y_min)
            && stan::length(alpha)))
        return 0.0;

      T_partials_return logp(0.0);

      check_not_nan(function, "Random variable", y);
      check_positive_finite(function, "Scale parameter", y_min);
      check_positive_finite(function, "Shape parameter", alpha);
      check_consistent_sizes(function,
                             "Random variable", y,
                             "Scale parameter", y_min,
                             "Shape parameter", alpha);

      if (!include_summand<propto, T_y, T_scale, T_shape>::value)
        return 0.0;

      VectorView<const T_y> y_vec(y);
      VectorView<const T_scale> y_min_vec(y_min);
      VectorView<const T_shape> alpha_vec(alpha);
      size_t N = max_size(y, y_min, alpha);

      for (size_t n = 0; n < N; n++) {
        if (y_vec[n] < y_min_vec[n])
          return LOG_ZERO;
      }

      OperandsAndPartials<T_y, T_scale, T_shape>
        operands_and_partials(y, y_min, alpha);

      VectorBuilder<include_summand<propto, T_y, T_shape>::value,
                    T_partials_return, T_y> log_y(length(y));
      if (include_summand<propto, T_y, T_shape>::value) {
        for (size_t n = 0; n < length(y); n++)
          log_y[n] = log(value_of(y_vec[n]));
      }

      VectorBuilder<contains_nonconstant_struct<T_y, T_shape>::value,
                    T_partials_return, T_y> inv_y(length(y));
      if (contains_nonconstant_struct<T_y, T_shape>::value) {
        for (size_t n = 0; n < length(y); n++)
          inv_y[n] = 1 / value_of(y_vec[n]);
      }

      VectorBuilder<include_summand<propto, T_scale, T_shape>::value,
                    T_partials_return, T_scale>
        log_y_min(length(y_min));
      if (include_summand<propto, T_scale, T_shape>::value) {
        for (size_t n = 0; n < length(y_min); n++)
          log_y_min[n] = log(value_of(y_min_vec[n]));
      }

      VectorBuilder<include_summand<propto, T_shape>::value,
                    T_partials_return, T_shape> log_alpha(length(alpha));
      if (include_summand<propto, T_shape>::value) {
        for (size_t n = 0; n < length(alpha); n++)
          log_alpha[n] = log(value_of(alpha_vec[n]));
      }

      for (size_t n = 0; n < N; n++) {
        const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
        if (include_summand<propto, T_shape>::value)
          logp += log_alpha[n];
        if (include_summand<propto, T_scale, T_shape>::value)
          logp += alpha_dbl * log_y_min[n];
        if (include_summand<propto, T_y, T_shape>::value)
          logp -= alpha_dbl * log_y[n] + log_y[n];

        if (!is_constant_struct<T_y>::value)
          operands_and_partials.d_x1[n] -= alpha_dbl * inv_y[n] + inv_y[n];
        if (!is_constant_struct<T_scale>::value)
          operands_and_partials.d_x2[n] += alpha_dbl / value_of(y_min_vec[n]);
        if (!is_constant_struct<T_shape>::value)
          operands_and_partials.d_x3[n]
            += 1 / alpha_dbl + log_y_min[n] - log_y[n];
      }
      return operands_and_partials.value(logp);
    }

    template <typename T_y, typename T_scale, typename T_shape>
    inline
    typename return_type<T_y, T_scale, T_shape>::type
    pareto_lpdf(const T_y& y, const T_scale& y_min, const T_shape& alpha) {
      return pareto_lpdf<false>(y, y_min, alpha);
    }

  }
}
#endif
