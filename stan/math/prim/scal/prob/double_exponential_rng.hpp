#ifndef STAN_MATH_PRIM_SCAL_PROB_DOUBLE_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_DOUBLE_EXPONENTIAL_RNG_HPP

#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/sign.hpp>

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom double exponential variate with the given location
     * and scale using the specified random number generator.
     *
     * @tparam RNG class of random number generator
     * @param mu location parameter
     * @param sigma positive scale parameter
     * @param rng random number generator
     * @return double exponential random variate
     * @throw std::domain_error if mu is infinite or sigma is nonpositive
     */
    template <class RNG>
    inline double
    double_exponential_rng(double mu,
                           double sigma,
                           RNG& rng) {
      static const char* function("double_exponential_rng");

      using boost::variate_generator;
      using boost::random::uniform_01;
      using std::log;
      using std::abs;

      check_finite(function, "Location parameter", mu);
      check_positive_finite(function, "Scale parameter", sigma);

      variate_generator<RNG&, uniform_01<> >
        rng_unit_01(rng, uniform_01<>());
      double a = 0;
      double laplaceRN = rng_unit_01();
      if (0.5 - laplaceRN > 0)
        a = 1.0;
      else if (0.5 - laplaceRN < 0)
        a = -1.0;
      return mu - sigma * a * log1m(2 * abs(0.5 - laplaceRN));
    }
  }
}
#endif
