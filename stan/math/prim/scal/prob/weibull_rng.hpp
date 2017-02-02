#ifndef STAN_MATH_PRIM_SCAL_PROB_WEIBULL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_WEIBULL_RNG_HPP

#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
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

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom Weibull variate with given shape and scale using
     * the specified random number generator.
     *
     * @tparam RNG class of random number generator
     * @param alpha shape parameter
     * @param sigma scale parameter
     * @param rng random number generator
     * @return Weibull random variate
     * @throw std::domain_error if alpha or sigma is nonpositive
     */
    template <class RNG>
    inline double
    weibull_rng(double alpha,
                double sigma,
                RNG& rng) {
      using boost::variate_generator;
      using boost::random::weibull_distribution;

      static const char* function("weibull_rng");

      check_positive_finite(function, "Shape parameter", alpha);
      check_positive_finite(function, "Scale parameter", sigma);

      variate_generator<RNG&, weibull_distribution<> >
        weibull_rng(rng, weibull_distribution<>(alpha, sigma));
      return weibull_rng();
    }

  }
}
#endif
