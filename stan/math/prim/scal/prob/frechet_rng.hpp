#ifndef STAN_MATH_PRIM_SCAL_PROB_FRECHET_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_FRECHET_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/VectorView.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    frechet_rng(double alpha,
                double sigma,
                RNG& rng) {
      using boost::variate_generator;
      using boost::random::weibull_distribution;

      static const char* function("frechet_rng");

      check_finite(function, "Shape parameter", alpha);
      check_positive(function, "Shape parameter", alpha);
      check_not_nan(function, "Scale parameter", sigma);
      check_positive(function, "Scale parameter", sigma);

      variate_generator<RNG&, weibull_distribution<> >
        weibull_rng(rng, weibull_distribution<>(alpha, 1.0/sigma));
      return 1.0 / weibull_rng();
    }

  }
}
#endif
