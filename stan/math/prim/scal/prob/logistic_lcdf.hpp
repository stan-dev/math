#ifndef STAN_MATH_PRIM_SCAL_PROB_LOGISTIC_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_LOGISTIC_LCDF_HPP

#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <stan/math/prim/scal/prob/logistic_log.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <cmath>
#include <limits>

namespace stan {
  namespace math {

    template <typename T_y, typename T_loc, typename T_scale>
    typename return_type<T_y, T_loc, T_scale>::type
    logistic_lcdf(const T_y& y, const T_loc& mu, const T_scale& sigma) {
      typedef typename stan::partials_return_type<T_y, T_loc, T_scale>::type
        T_partials_return;

      if (!(stan::length(y) && stan::length(mu) && stan::length(sigma)))
        return 0.0;

      static const char* function("logistic_lcdf");

      using boost::math::tools::promote_args;
      using std::log;
      using std::exp;

      T_partials_return P(0.0);

      check_not_nan(function, "Random variable", y);
      check_finite(function, "Location parameter", mu);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Random variable", y,
                             "Location parameter", mu,
                             "Scale parameter", sigma);

      scalar_seq_view<T_y> y_vec(y);
      scalar_seq_view<T_loc> mu_vec(mu);
      scalar_seq_view<T_scale> sigma_vec(sigma);
      size_t N = max_size(y, mu, sigma);

      operands_and_partials<T_y, T_loc, T_scale>
        ops_partials(y, mu, sigma);

      // Explicit return for extreme values
      // The gradients are technically ill-defined, but treated as zero
      for (size_t i = 0; i < stan::length(y); i++) {
        if (value_of(y_vec[i]) == -std::numeric_limits<double>::infinity())
          return ops_partials
            .build(-std::numeric_limits<double>::infinity());
      }

      for (size_t n = 0; n < N; n++) {
        // Explicit results for extreme values
        // The gradients are technically ill-defined, but treated as zero
        if (value_of(y_vec[n]) == std::numeric_limits<double>::infinity()) {
          continue;
        }

        const T_partials_return y_dbl = value_of(y_vec[n]);
        const T_partials_return mu_dbl = value_of(mu_vec[n]);
        const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
        const T_partials_return sigma_inv_vec = 1.0 / value_of(sigma_vec[n]);

        const T_partials_return Pn = 1.0 / (1.0 + exp(-(y_dbl - mu_dbl)
                                                      *sigma_inv_vec));
        P += log(Pn);

        if (!is_constant_struct<T_y>::value)
          ops_partials.edge1_.partials_[n]
            += exp(logistic_log(y_dbl, mu_dbl, sigma_dbl)) / Pn;
        if (!is_constant_struct<T_loc>::value)
          ops_partials.edge2_.partials_[n]
            += - exp(logistic_log(y_dbl, mu_dbl, sigma_dbl)) / Pn;
        if (!is_constant_struct<T_scale>::value)
          ops_partials.edge3_.partials_[n] += - (y_dbl - mu_dbl) * sigma_inv_vec
            * exp(logistic_log(y_dbl, mu_dbl, sigma_dbl)) / Pn;
      }
      return ops_partials.build(P);
    }

  }
}
#endif
