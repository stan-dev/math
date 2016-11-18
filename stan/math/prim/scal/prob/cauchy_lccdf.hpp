#ifndef STAN_MATH_PRIM_SCAL_PROB_CAUCHY_LCCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CAUCHY_LCCDF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/OperandsAndPartials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
  namespace math {

    template <typename T_y, typename T_loc, typename T_scale>
    typename return_type<T_y, T_loc, T_scale>::type
    cauchy_lccdf(const T_y& y, const T_loc& mu, const T_scale& sigma) {
      typedef typename stan::partials_return_type<T_y, T_loc, T_scale>::type
        T_partials_return;

      if ( !( stan::length(y) && stan::length(mu)
              && stan::length(sigma) ) )
        return 0.0;

      static const char* function("cauchy_lccdf");

      using boost::math::tools::promote_args;

      T_partials_return ccdf_log(0.0);

      check_not_nan(function, "Random variable", y);
      check_finite(function, "Location parameter", mu);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Random variable", y,
                             "Location parameter", mu,
                             "Scale Parameter", sigma);

      VectorView<const T_y> y_vec(y);
      VectorView<const T_loc> mu_vec(mu);
      VectorView<const T_scale> sigma_vec(sigma);
      size_t N = max_size(y, mu, sigma);

      OperandsAndPartials<T_y, T_loc, T_scale>
        operands_and_partials(y, mu, sigma);

      using std::atan;
      using std::log;

      for (size_t n = 0; n < N; n++) {
        const T_partials_return y_dbl = value_of(y_vec[n]);
        const T_partials_return mu_dbl = value_of(mu_vec[n]);
        const T_partials_return sigma_inv_dbl = 1.0 / value_of(sigma_vec[n]);
        const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
        const T_partials_return z = (y_dbl - mu_dbl) * sigma_inv_dbl;

        const T_partials_return Pn = 0.5 - atan(z) / pi();
        ccdf_log += log(Pn);

        const T_partials_return rep_deriv = 1.0 / (Pn * pi()
                                                   * (z * z * sigma_dbl
                                                      + sigma_dbl));
        if (!is_constant_struct<T_y>::value)
          operands_and_partials.d_x1[n] -= rep_deriv;
        if (!is_constant_struct<T_loc>::value)
          operands_and_partials.d_x2[n] += rep_deriv;
        if (!is_constant_struct<T_scale>::value)
          operands_and_partials.d_x3[n] += rep_deriv * z;
      }
      return operands_and_partials.value(ccdf_log);
    }

  }
}
#endif
