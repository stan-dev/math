#ifndef STAN_MATH_PRIM_SCAL_PROB_WEIBULL_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_WEIBULL_LCCDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/VectorView.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
  namespace math {

    template <typename T_y, typename T_shape, typename T_scale>
    typename return_type<T_y, T_shape, T_scale>::type
    weibull_lccdf(const T_y& y, const T_shape& alpha, const T_scale& sigma) {
      typedef typename stan::partials_return_type<T_y, T_shape, T_scale>::type
        T_partials_return;

      static const char* function("weibull_lccdf");

      using boost::math::tools::promote_args;
      using std::log;

      if (!(stan::length(y)
            && stan::length(alpha)
            && stan::length(sigma)))
        return 0.0;

      T_partials_return ccdf_log(0.0);
      check_nonnegative(function, "Random variable", y);
      check_positive_finite(function, "Shape parameter", alpha);
      check_positive_finite(function, "Scale parameter", sigma);

      OperandsAndPartials<T_y, T_shape, T_scale>
        operands_and_partials(y, alpha, sigma);

      VectorView<const T_y> y_vec(y);
      VectorView<const T_scale> sigma_vec(sigma);
      VectorView<const T_shape> alpha_vec(alpha);
      size_t N = max_size(y, sigma, alpha);
      for (size_t n = 0; n < N; n++) {
        const T_partials_return y_dbl = value_of(y_vec[n]);
        const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
        const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
        const T_partials_return pow_ = pow(y_dbl / sigma_dbl, alpha_dbl);

        ccdf_log -= pow_;

        if (!is_constant_struct<T_y>::value)
          operands_and_partials.d_x1[n] -= alpha_dbl / y_dbl * pow_;
        if (!is_constant_struct<T_shape>::value)
          operands_and_partials.d_x2[n] -= log(y_dbl / sigma_dbl) * pow_;
        if (!is_constant_struct<T_scale>::value)
          operands_and_partials.d_x3[n] += alpha_dbl / sigma_dbl * pow_;
      }
      return operands_and_partials.value(ccdf_log);
    }

  }
}
#endif
