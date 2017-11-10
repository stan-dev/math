#ifndef STAN_MATH_PRIM_SCAL_PROB_CAUCHY_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CAUCHY_RNG_HPP

#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom Cauchy variate for the given location and scale
     * using the specified random number generator.
     *
     * @tparam RNG type of random number generator
     * @param mu location parameter
     * @param sigma positive scale parameter
     * @param rng random number generator
     * @return Cauchy random variate
     * @throw std::domain_error if mu is infinite or sigma is nonpositive
     */
    template <class RNG>
    inline double
    cauchy_rng(double mu,
               double sigma,
               RNG& rng) {
      using boost::variate_generator;
      using boost::random::cauchy_distribution;

      static const char* function("cauchy_rng");

      check_finite(function, "Location parameter", mu);
      check_positive_finite(function, "Scale parameter", sigma);

      variate_generator<RNG&, cauchy_distribution<> >
        cauchy_rng(rng, cauchy_distribution<>(mu, sigma));
      return cauchy_rng();
    }

  }
}
#endif
