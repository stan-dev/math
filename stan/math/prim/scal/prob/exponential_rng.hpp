#ifndef STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    exponential_rng(double beta,
                    RNG& rng) {
      using boost::variate_generator;
      using boost::exponential_distribution;

      static const char* function("exponential_rng");

      check_positive_finite(function, "Inverse scale parameter", beta);

      variate_generator<RNG&, exponential_distribution<> >
        exp_rng(rng, exponential_distribution<>(beta));
      return exp_rng();
    }

  }
}
#endif
