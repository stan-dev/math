#ifndef STAN_MATH_PRIM_SCAL_PROB_GAMMA_LCDF_HPP
#define STAN_MATH_PRIM_SCAL_PROB_GAMMA_LCDF_HPP

#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>
#include <cmath>
#include <limits>

namespace stan {
  namespace math {

    template <typename T_y, typename T_shape, typename T_inv_scale>
    typename return_type<T_y, T_shape, T_inv_scale>::type
    gamma_lcdf(const T_y& y, const T_shape& alpha, const T_inv_scale& beta) {
      if (!(stan::length(y) && stan::length(alpha) && stan::length(beta)))
        return 0.0;
      typedef typename stan::partials_return_type<T_y, T_shape,
                                                  T_inv_scale>::type
        T_partials_return;

      static const char* function("gamma_lcdf");

      using boost::math::tools::promote_args;
      using std::exp;

      T_partials_return P(0.0);

      check_positive_finite(function, "Shape parameter", alpha);
      check_positive_finite(function, "Inverse scale parameter", beta);
      check_not_nan(function, "Random variable", y);
      check_nonnegative(function, "Random variable", y);
      check_consistent_sizes(function,
                             "Random variable", y,
                             "Shape parameter", alpha,
                             "Inverse scale parameter", beta);

      scalar_seq_view<T_y> y_vec(y);
      scalar_seq_view<T_shape> alpha_vec(alpha);
      scalar_seq_view<T_inv_scale> beta_vec(beta);
      size_t N = max_size(y, alpha, beta);

      operands_and_partials<T_y, T_shape, T_inv_scale>
        ops_partials(y, alpha, beta);

      // Explicit return for extreme values
      // The gradients are technically ill-defined, but treated as zero
      for (size_t i = 0; i < stan::length(y); i++) {
        if (value_of(y_vec[i]) == 0)
          return ops_partials.build(negative_infinity());
      }

      using boost::math::tgamma;
      using std::exp;
      using std::pow;
      using std::log;

      VectorBuilder<!is_constant_struct<T_shape>::value,
                    T_partials_return, T_shape> gamma_vec(stan::length(alpha));
      VectorBuilder<!is_constant_struct<T_shape>::value,
                    T_partials_return, T_shape>
        digamma_vec(stan::length(alpha));

      if (!is_constant_struct<T_shape>::value) {
        for (size_t i = 0; i < stan::length(alpha); i++) {
          const T_partials_return alpha_dbl = value_of(alpha_vec[i]);
          gamma_vec[i] = tgamma(alpha_dbl);
          digamma_vec[i] = digamma(alpha_dbl);
        }
      }

      for (size_t n = 0; n < N; n++) {
        // Explicit results for extreme values
        // The gradients are technically ill-defined, but treated as zero
        if (value_of(y_vec[n]) == std::numeric_limits<double>::infinity())
          return ops_partials.build(0.0);

        const T_partials_return y_dbl = value_of(y_vec[n]);
        const T_partials_return alpha_dbl = value_of(alpha_vec[n]);
        const T_partials_return beta_dbl = value_of(beta_vec[n]);

        const T_partials_return Pn = gamma_p(alpha_dbl, beta_dbl * y_dbl);

        P += log(Pn);

        if (!is_constant_struct<T_y>::value)
          ops_partials.edge1_.partials_[n] += beta_dbl * exp(-beta_dbl * y_dbl)
            * pow(beta_dbl * y_dbl, alpha_dbl - 1) / tgamma(alpha_dbl) / Pn;
        if (!is_constant_struct<T_shape>::value)
          ops_partials.edge2_.partials_[n]
            -= grad_reg_inc_gamma(alpha_dbl, beta_dbl
                                  * y_dbl, gamma_vec[n],
                                  digamma_vec[n]) / Pn;
        if (!is_constant_struct<T_inv_scale>::value)
          ops_partials.edge3_.partials_[n] += y_dbl * exp(-beta_dbl * y_dbl)
            * pow(beta_dbl * y_dbl, alpha_dbl - 1) / tgamma(alpha_dbl) / Pn;
      }
      return ops_partials.build(P);
    }

  }
}
#endif
